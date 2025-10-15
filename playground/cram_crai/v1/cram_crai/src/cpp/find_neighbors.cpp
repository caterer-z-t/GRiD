/*
COMPILATION COMMAND:
This compiles the program with optimization, includes necessary libraries for file I/O and compression
g++ \
-O3 \                                    # Optimization level 3 (fastest)
-Wall \                                  # Show all warnings
-std=c++03 \                            # Use C++03 standard
find_neighbors.cpp \                     # Source file
-I [path_to_eagle] \                    # Include directory for Eagle software
-I [path_to_boost] \                    # Include directory for Boost libraries
-Wl,-rpath,[boost_lib_path] \           # Runtime library path
-o find_neighbors \                      # Output executable name
-L [boost_lib_path] \                   # Library search path
-l boost_iostreams \                    # Link Boost iostreams (for compressed files)
-lz                                     # Link zlib (for compression)

USAGE:
./find_neighbors [batch_num] [total_batches] [max_z_range] [input_file] [output_prefix]

EXAMPLE:
./find_neighbors 0 1 2.0 input_zdepths.txt.gz output_prefix > log.txt

PURPOSE:
This program finds the nearest neighbors for genomic samples based on their depth coverage patterns.
It processes large datasets in batches to manage memory usage efficiently.
*/

// Standard C++ library includes
#include <iostream>    // For input/output operations (cout, cin, etc.)
#include <iomanip>     // For formatting output (precision, fixed point)
#include <sstream>     // For string stream operations
#include <fstream>     // For file operations
#include <vector>      // For dynamic arrays
#include <algorithm>   // For sorting and other algorithms
#include <cstdio>      // For C-style I/O functions (sprintf, sscanf)
#include <cmath>       // For mathematical functions

// Custom utility files (these must be in the same directory)
#include "FileUtils.cpp"   // Utilities for file handling (compressed files)
#include "StringUtils.cpp" // String manipulation utilities  
#include "Timer.cpp"       // Timer for performance monitoring

using namespace std;

// UTILITY FUNCTIONS

// Inline function to square a number (x^2)
// 'inline' suggests to compiler to replace function calls with actual code for speed
inline float sq(float x) { 
    return x*x; 
}

// Crop/clamp a z-score value to be within [-zMax, zMax] range
// This prevents extreme outliers from dominating distance calculations
float crop(float z, float zMax) {
    return min(zMax, max(-zMax, z));  // Clamp z to [-zMax, zMax]
}

// MAIN PROGRAM
int main(int argc, char *argv[]) {

    // Check if correct number of command line arguments provided
    if (argc != 6) {
        cerr << "ERROR: 5 arguments required" << endl;
        cerr << "- arg1: batch number (which batch to process)" << endl;
        cerr << "- arg2: total batches (how many batches total)" << endl;
        cerr << "- arg3: Max z range (maximum z-score value to allow)" << endl;
        cerr << "- arg4: ID_scale_zdepths file (input data file)" << endl;
        cerr << "- arg5: output prefix (prefix for output filename)" << endl;
        return 1;  // Exit with error code
    }
    
    // Parse command line arguments
    int b;                              // Current batch number (0-based)
    sscanf(argv[1], "%d", &b);         // Convert string to integer
    int B;                              // Total number of batches  
    sscanf(argv[2], "%d", &B);         // Convert string to integer
    const float fracR = 1;              // Fraction of regions to use (hardcoded to 1 = 100%)
    float zMax;                         // Maximum allowed z-score value
    sscanf(argv[3], "%f", &zMax);      // Convert string to float
    const char *dataFile = argv[4];     // Input file path
    const char *outPrefix = argv[5];    // Output file prefix

    // Print processing information
    cout << "Computing nearest neighbors for batch " << b << " mod " << B << endl;
    cout << "Cropping 'z-score' values to zMax = " << zMax << endl;

    Timer timer;  // Start timing the process

    // PHASE 1: READ FILE HEADER AND DETERMINE BATCH SIZE
    
    // Open compressed input file (handles .gz files automatically)
    FileUtils::AutoGzIfstream fin; 
    fin.openOrExit(dataFile);
    
    int N, R;           // N = number of individuals, R = number of genomic regions
    fin >> N >> R;      // Read first line: total individuals and regions
    
    // Calculate how many individuals will be in this specific batch
    int Nbatch = 0;
    for (int n = 0; n < N; n++)
        if (n % B == b)      // If individual n belongs to batch b
            Nbatch++;        // Count it
    cout << "the actual N_batch is " << Nbatch << endl;

    // Skip the header line containing "mu" (mean values)
    string line; 
    getline(fin, line);

    // PHASE 2: READ SIGMA-SQUARED RATIOS (VARIANCE INFORMATION)
    
    fin >> N >> R;  // Read dimensions again (file format requirement)
    vector<float> sigma2ratios(R);  // Store variance ratios for each region
    for (int r = 0; r < R; r++)
        fin >> sigma2ratios[r];     // Read variance ratio for region r

    // PHASE 3: FILTER REGIONS BASED ON VARIANCE
    // Remove regions with extreme variance (too high or too low)
    
    vector<float> sigma2sort = sigma2ratios;  // Copy for sorting
    sort(sigma2sort.begin(), sigma2sort.end()); // Sort from lowest to highest
    
    // Set thresholds for which regions to keep
    float sigma2min = sigma2sort[(int)(R*(1-fracR))]; // Minimum variance threshold
    float sigma2max = 1000;                           // Maximum variance threshold
    
    // Count how many regions we'll actually use
    int Ruse = 0, Rextreme = 0;
    for (int r = 0; r < R; r++) {
        if (sigma2ratios[r] >= sigma2min && sigma2ratios[r] <= sigma2max)
            Ruse++;           // Count usable regions
        if (sigma2ratios[r] > sigma2max)
            Rextreme++;       // Count extreme regions
    }
    
    // Report filtering results
    cout << "Removed " << Rextreme << " of " << R << " regions with sigma2ratio > " << sigma2max << endl;
    cout << "Keeping " << Ruse << " of " << R-Rextreme << " remaining regions with sigma2ratio >= " << sigma2min << endl;
    cout << "Reading data for " << Nbatch << " / " << N << " indivs in batch at " << R << " regions" << endl;

    // PHASE 4: ALLOCATE MEMORY FOR BATCH DATA
    
    vector<string> IDs(N);       // Store all individual IDs
    vector<float> scales(N);     // Store all individual scale factors
    float *zs = new float[R*(long long) Nbatch];  // Store z-scores for batch individuals only

    // PHASE 5: READ Z-SCORES FOR INDIVIDUALS IN THIS BATCH
    
    string ID; 
    float scale, z;
    for (int n = 0; n < N; n++) {
        fin >> ID >> scale;    // Read individual ID and scale factor
        IDs[n] = ID;           // Store ID for all individuals
        scales[n] = scale;     // Store scale for all individuals
        
        if (n % B == b) {      // If this individual belongs to our batch
            int i = n/B;       // Calculate local index within batch
            for (int r = 0; r < R; r++) {
                fin >> z;      // Read z-score for this region
                zs[r*Nbatch+i] = crop(z, zMax);  // Store cropped z-score
                // Array layout: regions are rows, batch individuals are columns
            }
        }
        else {
            getline(fin, line); // Skip z-scores for individuals not in this batch
        }
        
        if (n%100==0)
            cout << "." << flush;  // Progress indicator
    }
    fin.close();
    cout << endl << "Read data for " << Nbatch << " / " << N << " indivs in batch ("
         << timer.update_time() << " sec)" << endl;

    // PHASE 6: COMPUTE DISTANCES BETWEEN BATCH INDIVIDUALS AND ALL INDIVIDUALS
    
    // Allocate memory for distance matrix: [all individuals] x [batch individuals]
    float *dists = new float[N*(long long) Nbatch];
    memset(dists, 0, N*(long long) Nbatch*sizeof(dists[0])); // Initialize to zero
    
    // Re-read the file to compute distances
    fin.openOrExit(dataFile);
    getline(fin, line); // Skip first header line
    getline(fin, line); // Skip second header line
    
    // For each individual in the entire dataset
    for (int n = 0; n < N; n++) {
        fin >> ID >> scale; // Read ID and scale (already stored, just need to skip)
        
        // For each genomic region
        for (int r = 0; r < R; r++) {
            fin >> z;  // Read z-score for individual n at region r
            
            // Skip regions that were filtered out
            if (sigma2ratios[r] < sigma2min || sigma2ratios[r] > sigma2max) 
                continue;
            
            z = crop(z, zMax);  // Apply same cropping as batch data
            
            // Compute squared differences with all individuals in batch
            for (int i = 0; i < Nbatch; i++){
                dists[n*Nbatch+i] += sq(z - zs[r*Nbatch+i]);
                // Accumulate squared Euclidean distance
                // dists[n][i] = sum over regions of (z_n_r - z_i_r)^2
            }
        }
        
        if (n%100==0)
            cout << "." << flush;  // Progress indicator
    }
    fin.close();
    cout << endl << "Computed distances for " << Nbatch << " / " << N << " indivs in batch ("
         << timer.update_time() << " sec)" << endl;

    // PHASE 7: OUTPUT NEAREST NEIGHBORS
    
    // Create output filename
    char buf[200];
    sprintf(buf, "%s.zMax%g.txt.gz", outPrefix, zMax);
    
    // Open compressed output file
    FileUtils::AutoGzOfstream fout; 
    fout.openOrExit(buf);
    fout << std::setprecision(2) << std::fixed;  // Set output format
    
    // Vector to store (distance, individual_index) pairs for sorting
    vector<pair<float, int> > distIDs(N);
    
    // For each individual in the batch
    for (int i = 0; i < Nbatch; i++) {
        int n_i = i*B+b;  // Convert batch index back to global index
        fout << IDs[n_i] << "\t" << scales[n_i];  // Output ID and scale
        
        // Create pairs of (distance, individual_index) for sorting
        for (int n = 0; n < N; n++)
            distIDs[n] = make_pair(dists[n*Nbatch+i], n);
        
        // Set self-distance to infinity to exclude from neighbors
        distIDs[n_i].first = 1e9;
        
        // Sort by distance (ascending order)
        sort(distIDs.begin(), distIDs.end());

        // Output top 500 nearest neighbors
        const int N_OUTPUT = 500;
        for (int j = 0; j < N_OUTPUT; j++) {
            int n = distIDs[j].second;  // Get individual index
            // Output: neighbor_ID, neighbor_scale, normalized_distance
            fout << "\t" << IDs[n] << "\t" << scales[n] << "\t" << dists[n*Nbatch+i]/(2*Ruse);
            // Distance is normalized by 2*Ruse (twice the number of regions used)
        }
        fout << endl;
    }
    fout.close();
    cout << "Found neighbors and wrote output (" << timer.update_time() << " sec)" << endl;

    // CLEANUP: Free allocated memory
    delete[] dists;
    delete[] zs;
    
    return 0;  // Successful completion
}