/*
PROGRAM PURPOSE: Normalize sequencing depth data from mosdepth batch output files
- Reads depth data from multiple batch files
- Filters regions based on depth criteria and repeat overlaps
- Normalizes data both by individual and by region
- Outputs z-score normalized depth values

COMPILATION:
g++ -O3 -Wall \
    -std=c++03 normalize_mosdepth.cpp \
    -I$SOFTWARE/Eagle_v2.4.1/src/ \
    -I$SOFTWARE/boost/include \
    -L$SOFTWARE/boost/lib \
    -lboost_iostreams -lz -o normalize_mosdepth

USAGE:
./normalize_mosdepth \
<prefix> \                    # prefix of mosdepth batch files
<repeat_mask_bed> \          # bed file with repeat regions to exclude
<example_regions_file> \     # example .regions.bed.gz file for region coordinates
<total_sample_size> \        # total number of samples
<output_path>                # output file path (.txt.gz)
*/

#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cstdio>
#include <cmath>
#include <set>
#include "FileUtils.cpp"    // Custom file utilities (handles gzipped files)
#include "StringUtils.cpp"  // Custom string utilities
#include "Timer.cpp"        // Custom timer for performance monitoring

using namespace std;

int main(int argc, char *argv[]) {

  // CHECK COMMAND LINE ARGUMENTS
  if (argc != 6) {
    cerr << "ERROR: 5 arguments required" << endl;
    cerr << "- arg1: prefix of mosdepth input (no more than 170 characters), used in <prefix>_batch_<batchnumber>.txt.gz" << endl;
    cerr << "- arg2: bed file path e.g. /path/to/repeat_mask_list.hg38.ucsc_bed" << endl;
    cerr << "- arg3: example input e.g. /path/to/name_regions.bed.gz" << endl;
    cerr << "- arg4: N_sample(int)" << endl;
    cerr << "- arg5: output path e.g. /path/to/ID_scale_zdepths.txt.gz" << endl;
    return 1;
  }
  
  // PARSE COMMAND LINE ARGUMENTS
  const char *mosdepth_prefix = argv[1];  // Prefix for batch files
  const char *bed_source = argv[2];       // Repeat mask bed file
  const char *example_output = argv[3];   // Example regions file for coordinates
  int MAX_N; sscanf(argv[4], "%d", &MAX_N); // Maximum number of samples
  const char *output_path = argv[5];      // Output file path

  Timer timer; // Start timing execution
  FileUtils::AutoGzIfstream fin; // File input stream that handles gzipped files

  /**************************************************************************
   * PHASE 1: DECIDE WHICH REGIONS TO EXTRACT
   * Read first 10 batches (~250 individuals) to calculate mean depths
   * and select regions with reasonable depth (20-100x coverage)
   **************************************************************************/
  
  long long R = 0; // Total number of regions
  vector <double> mean_depths; // Running sum of depths for each region
  int N = 0; // Number of individuals processed so far
  
  // Process first 10 batches to get representative depth data
  for (int batch = 0; batch < 10; batch++) {
    char buf[200];
    // Create filename: <prefix>_batch_<batchnumber>.txt.gz
    sprintf(buf, "%s_batch_%d.txt.gz", mosdepth_prefix, batch+1);
    fin.openOrExit(buf);
    
    string ID; // Individual/sample identifier (not used in this phase)
    
    // For the first individual, initialize the regions count and mean_depths vector
    if (R == 0) {
      N++; // Count this individual
      string line;
      getline(fin, line); // Read entire first line
      istringstream iss(line); // Parse the line
      iss >> ID; // Extract sample ID
      
      int depth;
      // Count regions by reading all depth values from first line
      while (iss >> depth) {
        mean_depths.push_back(depth*0.01); // Convert to actual depth (mosdepth outputs 100x)
        R++; // Count total regions
      }
    }
    
    // Process remaining individuals in this batch
    while (fin >> ID) {
      N++; // Count this individual
      int depth;
      // Read depth for each region and add to running sum
      for (int r = 0; r < R; r++) {
        fin >> depth;
        mean_depths[r] += depth*0.01; // Add to running sum, convert from 100x
      }
    }
    fin.close();
    cout << "Read batch " << batch << endl;
  }

  // FILTER REGIONS BASED ON MEAN DEPTH
  vector <bool> extract(R); // Which regions to include in analysis
  long long Rextract = 0; // Number of regions that pass depth filter
  
  for (int r = 0; r < R; r++) {
    mean_depths[r] /= N; // Convert sum to mean
    // Keep regions with reasonable coverage (20-100x)
    if (20 <= mean_depths[r] && mean_depths[r] <= 100) {
      extract[r] = true;
      Rextract++;
    }
    // extract[r] remains false (default) for regions outside this range
  }

  cout << "Read " << N << " indivs (" << timer.update_time() << " sec)" << endl;
  cout << "Extracting " << Rextract << " / " << R << " regions" << endl;

  /**************************************************************************
   * PHASE 2: EXCLUDE REGIONS OVERLAPPING REPEATS/VNTRs
   * Read repeat mask file and mark 1kb regions that overlap repeats
   **************************************************************************/

  // Create a 2D boolean array: overlapsVNTR[chromosome][kilobase_position]
  // This tracks which 1kb windows overlap with repeat elements
  vector < vector <bool> > overlapsVNTR(23, vector <bool> (300000));
  fin.openOrExit(bed_source);
  
  int numVNTRs = 0, totLen = 0; // Statistics counters
  string chrStr;
  
  // Define which chromosomes to process (currently only chromosome 6)
  std::set<int> valid_chr;
  std::set<int>::iterator it = valid_chr.begin();
  it = valid_chr.insert(it, 6); // Only process chromosome 6
  
  // Uncomment these lines to process all autosomes (1-22):
  // for (int i = 1; i < 23; i++) {
  //     it = valid_chr.insert(it, i);
  // }

  // Read repeat regions from BED file
  while (fin >> chrStr) {
    if (chrStr == "chrX") continue; // Skip X chromosome
    
    int chr, bpStart, bpEnd, bpLen; 
    string vntrStr;
    
    // Parse chromosome number from string like "chr6"
    sscanf(chrStr.c_str(), "chr%d", &chr);
    
    // Skip chromosomes we're not interested in
    if (valid_chr.count(chr) == 0) continue;
    
    numVNTRs++;
    // Read the rest of the BED line: start, end, name, length
    fin >> bpStart >> bpEnd >> vntrStr >> bpLen;
    
    // Mark all 1kb windows that overlap this repeat region
    for (int kb = bpStart/1000; kb <= bpEnd/1000; kb++)
      overlapsVNTR[chr][kb] = true;
    
    totLen += bpEnd - bpStart; // Track total repeat length
  }
  fin.close();
  cout << "Read " << numVNTRs << " autosomal VNTRs spanning " << totLen*1e-6 << " Mb" << endl;

  // UPDATE EXTRACT LIST TO EXCLUDE REPEAT-OVERLAPPING REGIONS
  int Roverlap = 0, RextractOverlap = 0; // Statistics counters
  fin.openOrExit(example_output);
  
  // Read region coordinates and check for repeat overlap
  for (int r = 0; r < R; r++) {
    int chr, bpStart, bpEnd; 
    double dep; // Depth (not used here, just read to advance file pointer)
    fin >> chr >> bpStart >> bpEnd >> dep;
    
    // Skip regions on chromosomes we're not processing
    if (valid_chr.count(chr) == 0) {
      r--; // Don't count this region, repeat this iteration
      continue;
    }
    
    // Check if this region overlaps a repeat (check the 1kb window it falls in)
    if (overlapsVNTR[chr][bpStart/1000]) {
      Roverlap++; // Count total overlapping regions
      if (extract[r]) { // If this region was previously selected for extraction
        RextractOverlap++; // Count overlapping regions in extract set
        extract[r] = false; // Exclude it now
        Rextract--; // Decrease count of regions to extract
      }
    }
    
    // Print info about the last region for verification
    if (r==R-1) cout << "Last region: " << chr << ":" << bpStart << "-" << bpEnd << endl;
  }
  fin.close();
  cout << "Excluding " << Roverlap << " / " << R << " regions overlapping VNTRs" << endl;
  cout << "Excluded " << RextractOverlap << " in extract set; " << Rextract << " left" << endl;

  /**************************************************************************
   * PHASE 3: EXTRACT DATA FOR ALL INDIVIDUALS
   * Read all batch files and normalize each individual's data
   **************************************************************************/

  N = 0; // Reset counter for all individuals
  const double batch_size = 25.0; // Assumed batch size for calculating number of batches
  const int batch_num = ceil(MAX_N/batch_size); // Total number of batch files
  
  // Allocate memory for depth data: MAX_N individuals × Rextract regions
  float *depths = new float[MAX_N*Rextract];
  vector <string> IDs(MAX_N); // Store individual IDs
  vector <float> scales(MAX_N); // Store scaling factors for each individual

  // Process all batch files
  for (int batch = 0; batch < batch_num; batch++) {
    // Skip certain problematic batches (commented out - adjust as needed)
    //if (376<=batch && batch<=423 || 1000<=batch && batch<1700) continue;
    
    char buf[200];
    sprintf(buf, "%s_batch_%d.txt.gz", mosdepth_prefix, batch+1);
    fin.openOrExit(buf);
    
    string ID;
    // Process each individual in this batch
    while (fin >> ID) {
      // Get pointer to this individual's row in the depths array
      float *depthsRow = depths + N*Rextract;
      IDs[N] = ID; // Store the individual ID
      
      int rSub = 0; // Index for extracted regions only
      float mean_depthExtract = 0; // Running sum for mean calculation
      
      int depth;
      // Read depth for each region
      for (int r = 0; r < R; r++) {
        fin >> depth;
        // Check for missing data
        if (isnan(depth)){
          cout << "na: depth: " << r << batch << endl;
        }
        
        // If this region passed our filters, store it
        if (extract[r]) {
          depthsRow[rSub++] = depth; // Store depth value
          mean_depthExtract += depth; // Add to sum for mean
        }
      }
      
      // NORMALIZE THIS INDIVIDUAL'S DATA
      mean_depthExtract /= Rextract; // Calculate mean depth for this individual
      if (isnan(mean_depthExtract)){
        cout << "na: mean_depthExtract: " << batch << endl;
      }
      
      scales[N] = mean_depthExtract; // Store the scaling factor
      float invScale = 1 / mean_depthExtract; // Calculate inverse for normalization
      if (isnan(invScale)){
        cout << "na: invScale: "<< batch <<endl;
      }
      
      // Normalize each depth value by dividing by individual's mean depth
      for (int s = 0; s < Rextract; s++){
        depthsRow[s] *= invScale; // Now each value is relative to individual's mean
        if (isnan(depthsRow[s])){
          cout << "na: depthsRow[s]: " << s << batch <<endl;
        }
      }
      N++; // Count this individual
    }
    fin.close();
    cout << "Read batch " << batch << " (" << timer.update_time() << " sec)" << endl;
  }
  cout << "Read " << N << " indivs; normalizing by region" << endl;

  /**************************************************************************
   * PHASE 4: NORMALIZE EACH REGION ACROSS ALL INDIVIDUALS
   * Calculate z-scores for each region
   **************************************************************************/

  float x[N]; // Temporary array to hold values for one region across all individuals
  const float ratioMult = 100; // Multiplier for more readable output
  vector <float> mus(Rextract), sigma2s(Rextract), sigma2ratios(Rextract);

  // Process each region
  for (int s = 0; s < Rextract; s++) {
    // Extract values for this region across all individuals
    for (int n = 0; n < N; n++) {
      x[n] = depths[n*Rextract+s];
      if(isnan(x[n])){
        cout << "current info x: NA: n=" << n  << " s="<< s << endl;
      }
    }
    
    // CALCULATE MEAN AND VARIANCE FOR THIS REGION
    float mu = 0, s2 = 0;
    
    // Calculate mean
    for (int n = 0; n < N; n++) mu += x[n];
    mu /= N;
    
    // Calculate variance
    for (int n = 0; n < N; n++) s2 += (x[n]-mu)*(x[n]-mu);
    s2 /= (N-1); // Sample variance (N-1 denominator)
    
    // Store statistics for this region
    mus[s] = mu; 
    sigma2s[s] = s2; 
    sigma2ratios[s] = ratioMult*s2/mu; // Variance-to-mean ratio (scaled)
    
    // NORMALIZE VALUES FOR THIS REGION
    // Using (x-mu)/sqrt(mu) instead of standard z-score (x-mu)/sigma
    // This is a variance-stabilizing transformation
    float invRootMean = 1 / sqrtf(mu);
    for (int n = 0; n < N; n++) {
      depths[n*Rextract+s] = (x[n]-mu) * invRootMean;
    }
  }
  cout << "Normalized by region (" << timer.update_time() << " sec)" << endl;

  /**************************************************************************
   * PHASE 5: SELECT HIGH-VARIANCE REGIONS
   * Keep only regions with high variance-to-mean ratios (top 10%)
   **************************************************************************/

  // Sort variance-to-mean ratios to find thresholds
  sort(sigma2ratios.begin(), sigma2ratios.end());
  float sigma2ratioMin = sigma2ratios[(int) (0.9*Rextract)]; // 90th percentile threshold
  
  // Select regions above threshold
  vector <bool> want(Rextract); 
  int Rwant = 0;
  for (int s = 0; s < Rextract; s++) {
    if (ratioMult*sigma2s[s]/mus[s] > sigma2ratioMin) {
      want[s] = true;
      Rwant++;
    }
  }
  
  cout << "Restricting to " << Rwant << " regions with sigma2ratio > " << sigma2ratioMin << endl;
  
  // Calculate scaling factor based on median variance ratio
  float sigma2ratioMedian = sigma2ratios[Rextract/2];
  cout << "Rescaling to approximate z-scores based on median sigma2ratio = " << sigma2ratioMedian << endl;
  float scale = 1/sqrtf(sigma2ratioMedian/ratioMult);

  /**************************************************************************
   * PHASE 6: WRITE OUTPUT
   * Output format:
   * Line 1: N_samples N_regions mean1 mean2 ... (means for each region)
   * Line 2: N_samples N_regions var_ratio1 var_ratio2 ... (variance ratios)
   * Line 3+: sample_ID scale depth1 depth2 ... (normalized depths for each sample)
   **************************************************************************/

  FileUtils::AutoGzOfstream fout; 
  fout.openOrExit(output_path);
  fout << std::setprecision(3) << std::fixed; // Set output precision

  // WRITE HEADER LINE 1: Sample count, region count, and region means
  fout << N << "\t" << Rwant;
  for (int s = 0; s < Rextract; s++)
    if (want[s]) // Only output selected regions
      fout << "\t" << mus[s];
  fout << endl;

  // WRITE HEADER LINE 2: Sample count, region count, and variance ratios
  fout << N << "\t" << Rwant;
  for (int s = 0; s < Rextract; s++)
    if (want[s])
      fout << "\t" << ratioMult*sigma2s[s]/mus[s];
  fout << endl;

  // WRITE DATA LINES: Individual data
  fout << std::setprecision(2) << std::fixed; // Reduce precision for data
  for (int n = 0; n < N; n++) {
    // Output: ID, individual scaling factor, normalized depths
    fout << IDs[n] << "\t" << 0.01f*scales[n]; // Convert scale back to original units
    for (int s = 0; s < Rextract; s++)
      if (want[s])
        fout << "\t" << scale*depths[n*Rextract+s]; // Apply final scaling
    fout << endl;
  }
  
  fout.close();
  cout << "Wrote output (" << timer.update_time() << " sec)" << endl;

  // CLEANUP
  delete[] depths; // Free allocated memory
  
  return 0;
}