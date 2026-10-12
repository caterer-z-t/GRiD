// Compilation:
// g++ -O3 -Wall -std=c++14 \
//     normalize_mosdepth.cpp \
//     -I $SOFTWARE/Eagle_v2.4.1/src \
//     -I $SOFTWARE/boost/install/include \
//     -L $SOFTWARE/boost/install/lib \
//     -Wl,-rpath,$SOFTWARE/boost/install/lib \
//     -o ./normalize_mosdepth \
//     -lboost_iostreams -lz

#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <set>

#include "FileUtils.cpp"
#include "StringUtils.cpp"
#include "Timer.cpp"

using namespace std;

int main(int argc, char* argv[]) {
  if (argc < 6) {
    cerr << "Usage: " << argv[0] 
         << " <combined_bins_gz> <bed_mask_file> <example_regions_bed_gz> <num_samples> <output_gz> <min_depth: default=20> <max_depth: default=100>" 
         << endl;
    return 1;
  }

  string combined_bins_gz = argv[1];
  string bed_mask_file = argv[2];
  string example_regions_bed_gz = argv[3];
  int MAX_N = atoi(argv[4]);
  double min_depth = (argc > 6) ? atof(argv[6]) : 20.0;
  double max_depth = (argc > 7) ? atof(argv[7]) : 100.0;
  string output_gz = argv[5];

  Timer timer;
  FileUtils::AutoGzIfstream fin;

  /***** 1. DECIDE WHICH REGIONS TO EXTRACT *****/
  long long R = 0; // number of regions
  vector<double> mean_depths;
  int N = 0;

  fin.openOrExit(combined_bins_gz.c_str());

  string line;
  if (getline(fin, line)) {
    N++;
    istringstream iss(line);
    string sample_id;
    iss >> sample_id;
    int depth;
    while (iss >> depth) {
      mean_depths.push_back(depth * 0.01);
      R++;
    }
  }

  if (R == 0) {
    cerr << "ERROR: Failed to read regions from first line of " << combined_bins_gz << endl;
    return 1;
  }

  while (getline(fin, line)) {
    if (line.empty()) continue;
    N++;
    istringstream iss(line);
    string sample_id;
    iss >> sample_id;
    int depth;
    for (long long r = 0; r < R; r++) {
      if (iss >> depth) {
        mean_depths[r] += depth * 0.01;
      }
    }
  }
  fin.close();
  cout << "Read combined file with " << N << " samples across " << R << " regions." << endl;

  vector<bool> extract(R, false);
  long long Rextract = 0;
  for (long long r = 0; r < R; r++) {
    mean_depths[r] /= N;
    if (min_depth <= mean_depths[r] && mean_depths[r] <= max_depth) {
      extract[r] = true;
      Rextract++;
    }
  }

  cout << "Read " << N << " indivs (" << timer.update_time() << " sec)" << endl;
  cout << "Extracting " << Rextract << " / " << R << " regions based on mean depth [20, 100]" << endl;

  /***** 2. EXCLUDE REGIONS OVERLAPPING VNTRs *****/
  vector<vector<bool> > overlapsVNTR(23, vector<bool>(300000, false));
  fin.openOrExit(bed_mask_file.c_str());
  int numVNTRs = 0;
  long long totLen = 0;
  string chrStr;

  set<int> valid_chr;
  for (int c = 1; c <= 22; c++) valid_chr.insert(c);

  while (fin >> chrStr) {
    if (chrStr == "chrX" || chrStr == "X" || chrStr == "chrY" || chrStr == "Y") continue;
    numVNTRs++;
    int chr = 0, bpStart = 0, bpEnd = 0, bpLen = 0;
    string vntrStr;
    
    if (chrStr.find("chr") == 0) {
      sscanf(chrStr.c_str(), "chr%d", &chr);
    } else {
      sscanf(chrStr.c_str(), "%d", &chr);
    }

    fin >> bpStart >> bpEnd >> vntrStr >> bpLen;
    if (valid_chr.count(chr)) {
      for (int kb = bpStart / 1000; kb <= bpEnd / 1000; kb++) {
        if (kb < 300000) overlapsVNTR[chr][kb] = true;
      }
      totLen += (bpEnd - bpStart);
    }
  }
  fin.close();
  cout << "Read " << numVNTRs << " autosomal VNTRs spanning " << totLen * 1e-6 << " Mb" << endl;

  // Exclude overlapping regions
  int Roverlap = 0, RextractOverlap = 0;
  fin.openOrExit(example_regions_bed_gz.c_str());
  for (long long r = 0; r < R; r++) {
    string chr_raw;
    int bpStart, bpEnd;
    double dep;
    if (!(fin >> chr_raw >> bpStart >> bpEnd >> dep)) break;

    int chr = 0;
    if (chr_raw.find("chr") == 0) {
      sscanf(chr_raw.c_str(), "chr%d", &chr);
    } else {
      sscanf(chr_raw.c_str(), "%d", &chr);
    }

    if (valid_chr.count(chr) && overlapsVNTR[chr][bpStart / 1000]) {
      Roverlap++;
      if (extract[r]) {
        RextractOverlap++;
        extract[r] = false;
        Rextract--;
      }
    }
  }
  fin.close();
  cout << "Excluding " << Roverlap << " / " << R << " regions overlapping VNTRs" << endl;
  cout << "Excluded " << RextractOverlap << " in extract set; " << Rextract << " left" << endl;

  if (Rextract == 0) {
    cerr << "ERROR: 0 regions remaining to extract! Check bed mask or input thresholds." << endl;
    return 1;
  }

  /***** 3. EXTRACT DATA & NORMALIZE BY INDIVIDUAL *****/
  float* depths = new float[(long long)N * Rextract];
  vector<string> IDs(N);
  vector<float> scales(N);

  int current_n = 0;
  fin.openOrExit(combined_bins_gz.c_str());

  while (getline(fin, line) && current_n < N) {
    if (line.empty()) continue;
    istringstream iss(line);
    string sample_id;
    iss >> sample_id;
    IDs[current_n] = sample_id;

    float* depthsRow = depths + (long long)current_n * Rextract;
    long long rSub = 0;
    float mean_depthExtract = 0;
    int depth;

    for (long long r = 0; r < R; r++) {
      if (iss >> depth) {
        if (extract[r] && rSub < Rextract) {
          depthsRow[rSub++] = static_cast<float>(depth);
          mean_depthExtract += depth;
        }
      }
    }

    if (rSub > 0) {
      mean_depthExtract /= rSub;
    } else {
      mean_depthExtract = 1.0f;
    }

    scales[current_n] = mean_depthExtract;
    float invScale = (mean_depthExtract > 0.0f) ? (1.0f / mean_depthExtract) : 1.0f;

    for (long long s = 0; s < Rextract; s++) {
      depthsRow[s] *= invScale;
    }
    current_n++;
  }
  fin.close();
  cout << "Read " << current_n << " indivs; normalizing by region" << endl;

  /***** 4. NORMALIZE BY REGION ACROSS INDIVIDUALS *****/
  vector<float> x(N);
  const float ratioMult = 100.0f;
  vector<float> mus(Rextract), sigma2s(Rextract), sigma2ratios(Rextract);

  for (long long s = 0; s < Rextract; s++) {
    for (int n = 0; n < N; n++) x[n] = depths[(long long)n * Rextract + s];
    float mu = 0, s2 = 0;
    for (int n = 0; n < N; n++) mu += x[n];
    mu /= N;
    for (int n = 0; n < N; n++) s2 += (x[n] - mu) * (x[n] - mu);
    s2 /= (N > 1 ? (N - 1) : 1);

    mus[s] = mu;
    sigma2s[s] = s2;
    sigma2ratios[s] = (mu > 0) ? (ratioMult * s2 / mu) : 0;

    float invRootMean = (mu > 0) ? (1.0f / sqrtf(mu)) : 1.0f;
    for (int n = 0; n < N; n++) {
      depths[(long long)n * Rextract + s] = (x[n] - mu) * invRootMean;
    }
  }
  cout << "Normalized by region (" << timer.update_time() << " sec)" << endl;

  vector<float> sorted_ratios = sigma2ratios;
  sort(sorted_ratios.begin(), sorted_ratios.end());
  float sigma2ratioMin = sorted_ratios[(int)(0.9 * Rextract)];
  
  vector<bool> want(Rextract, false);
  int Rwant = 0;
  for (long long s = 0; s < Rextract; s++) {
    if (mus[s] > 0 && (ratioMult * sigma2s[s] / mus[s]) > sigma2ratioMin) {
      want[s] = true;
      Rwant++;
    }
  }

  cout << "Restricting to " << Rwant << " regions with sigma2ratio > " << sigma2ratioMin << endl;
  float sigma2ratioMedian = sorted_ratios[Rextract / 2];
  cout << "Rescaling to approximate z-scores based on median sigma2ratio = " << sigma2ratioMedian << endl;
  float scale = (sigma2ratioMedian > 0) ? (1.0f / sqrtf(sigma2ratioMedian / ratioMult)) : 1.0f;

  /***** 5. WRITE OUTPUT *****/
  FileUtils::AutoGzOfstream fout;
  fout.openOrExit(output_gz.c_str());
  fout << std::setprecision(3) << std::fixed;

  fout << N << "\t" << Rwant;
  for (long long s = 0; s < Rextract; s++) if (want[s]) fout << "\t" << mus[s];
  fout << endl;

  fout << N << "\t" << Rwant;
  for (long long s = 0; s < Rextract; s++) if (want[s]) fout << "\t" << ratioMult * sigma2s[s] / mus[s];
  fout << endl;

  fout << std::setprecision(2) << std::fixed;
  for (int n = 0; n < N; n++) {
    fout << IDs[n] << "\t" << 0.01f * scales[n];
    for (long long s = 0; s < Rextract; s++) {
      if (want[s]) fout << "\t" << scale * depths[(long long)n * Rextract + s];
    }
    fout << endl;
  }
  fout.close();
  cout << "Wrote output (" << timer.update_time() << " sec)" << endl;

  delete[] depths;
  return 0;
}