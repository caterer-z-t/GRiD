#!/usr/bin/env python3
"""
Python conversion of normalize_mosdepth_inflow_rewritten.cpp
Normalizes individual mosdepth .regions.bed.gz files for VNTR analysis

Usage:
python normalize_mosdepth.py <mosdepth_dir> <repeat_mask_bed> <n_samples> <output_file>
"""

import sys
import gzip
import math
import numpy as np
import time
from collections import defaultdict
import os

# Get the file path of the current script
current_file_path = os.path.abspath(__file__)
current_dir = os.path.dirname(current_file_path)
utils_path = os.path.join(current_dir, "utils", "utils.py")

# Import all contents from utils/utils.py
sys.path.insert(0, os.path.join(current_dir, "utils"))

from utils.utils import *


def determine_regions_to_extract(mosdepth_prefix, target_chr=6, max_batches=10):
    """
    Analyze first 10 batches to determine which regions to extract
    based on mean depth criteria (20-100x coverage after scaling by 0.01)
    """
    log_message("Analyzing batches to determine regions to extract...")
    
    R = 0  # number of regions
    mean_depths = []
    N = 0  # number of individuals processed
    
    for batch in range(max_batches):
        batch_file = f"{mosdepth_prefix}_batch_{batch + 1}.txt.gz"
        
        if not os.path.exists(batch_file):
            log_message(f"Warning: Batch file {batch_file} not found, continuing with available batches")
            continue
        
        try:
            with read_gzipped_file(batch_file) as f:
                for line in f:
                    if not line.strip():
                        continue
                    
                    parts = line.strip().split('\t')
                    if len(parts) < 2:
                        continue
                    
                    sample_id = parts[0]
                    depths = [int(d) for d in parts[1:]]
                    
                    # Initialize on first sample
                    if R == 0:
                        R = len(depths)
                        mean_depths = [0.0] * R
                    
                    if len(depths) != R:
                        log_message(f"Warning: Sample {sample_id} has {len(depths)} depths, expected {R}")
                        continue
                    
                    # Add depths (convert by 0.01 like C++ version)
                    for r, depth in enumerate(depths):
                        mean_depths[r] += depth * 0.01
                    
                    N += 1
        
        except Exception as e:
            log_message(f"Error reading batch {batch}: {e}")
            continue
        
        log_message(f"Read batch {batch}")
    
    # Calculate mean depths and determine extraction criteria
    extract_regions = []
    num_extract = 0
    
    for r in range(R):
        mean_depth = mean_depths[r] / N if N > 0 else 0
        # Use criteria: 20 <= mean_depth <= 100 (after 0.01 scaling)
        if 20 <= mean_depth <= 100:
            extract_regions.append(True)
            num_extract += 1
        else:
            extract_regions.append(False)
    
    log_message(f"Read {N} individuals")
    log_message(f"Extracting {num_extract} / {R} regions")
    
    return extract_regions, mean_depths, N

def load_repeat_mask(bed_file, target_chr=6):
    """Load repeat mask regions for specified chromosome"""
    log_message(f"Loading repeat mask from {bed_file}")
    
    overlaps_vntr = defaultdict(lambda: defaultdict(bool))
    num_vntrs = 0
    total_length = 0
    
    try:
        with open(bed_file, 'r') as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                
                parts = line.strip().split('\t')
                if len(parts) < 3:
                    continue
                
                chr_num = parse_chromosome(parts[0])
                if chr_num != target_chr:
                    continue
                
                start = int(parts[1])
                end = int(parts[2])
                
                num_vntrs += 1
                length = end - start
                total_length += length
                
                # Mark 1kb regions that overlap
                for kb in range(start // 1000, (end // 1000) + 1):
                    overlaps_vntr[chr_num][kb] = True
        
        log_message(f"Read {num_vntrs} autosomal VNTRs spanning {total_length/1e6:.2f} Mb")
        return overlaps_vntr
    
    except FileNotFoundError:
        log_message(f"Warning: Repeat mask file {bed_file} not found")
        return defaultdict(lambda: defaultdict(bool))

def update_extraction_list(example_regions_file, extract_regions, overlaps_vntr, target_chr=6):
    """Update extraction list to exclude regions overlapping VNTRs"""
    log_message("Updating extraction list to exclude VNTR overlaps...")
    
    regions_overlap = 0
    extract_overlap = 0
    
    try:
        with read_gzipped_file(example_regions_file) as f:
            r = 0
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                
                parts = line.strip().split('\t')
                if len(parts) < 4:
                    continue
                
                chr_num = parse_chromosome(parts[0])
                start = int(parts[1])
                end = int(parts[2])
                
                # Skip non-target chromosomes but don't increment counter
                if chr_num != target_chr:
                    continue
                
                if r >= len(extract_regions):
                    break
                
                # Check if region overlaps VNTR
                if overlaps_vntr[chr_num][start // 1000]:
                    regions_overlap += 1
                    if extract_regions[r]:
                        extract_overlap += 1
                        extract_regions[r] = False
                
                r += 1
        
        num_extract = sum(extract_regions)
        log_message(f"Excluding {regions_overlap} / {len(extract_regions)} regions overlapping VNTRs")
        log_message(f"Excluded {extract_overlap} in extract set; {num_extract} left")
    
    except Exception as e:
        log_message(f"Error processing example regions file: {e}")
    
    return extract_regions

def load_all_samples(mosdepth_prefix, extract_regions, max_n_samples):
    """Load depth data for all samples from batch files"""
    log_message("Loading all sample data...")
    
    num_extract = sum(extract_regions)
    if num_extract == 0:
        log_message("Error: No regions selected for extraction!")
        return [], [], np.array([])
    
    batch_size = 25.0
    batch_num = math.ceil(max_n_samples / batch_size)
    
    all_ids = []
    all_scales = []
    all_depths = []
    
    for batch in range(batch_num):
        batch_file = f"{mosdepth_prefix}_batch_{batch + 1}.txt.gz"
        
        if not os.path.exists(batch_file):
            log_message(f"Warning: Batch file {batch_file} not found, skipping")
            continue
        
        try:
            with read_gzipped_file(batch_file) as f:
                for line in f:
                    if not line.strip():
                        continue
                    
                    parts = line.strip().split('\t')
                    if len(parts) < 2:
                        continue
                    
                    sample_id = parts[0]
                    depths = [int(d) for d in parts[1:]]
                    
                    if len(depths) != len(extract_regions):
                        log_message(f"Warning: Sample {sample_id} has {len(depths)} depths, expected {len(extract_regions)}")
                        continue
                    
                    # Extract depths for selected regions
                    extract_depths = []
                    for r, depth in enumerate(depths):
                        if extract_regions[r]:
                            extract_depths.append(float(depth))
                    
                    if not extract_depths:
                        continue
                    
                    # Calculate mean depth and scale factor
                    mean_depth = np.mean(extract_depths)
                    scale_factor = mean_depth
                    
                    # Normalize by individual mean
                    if mean_depth > 0:
                        normalized_depths = np.array(extract_depths) / mean_depth
                    else:
                        normalized_depths = np.array(extract_depths)
                    
                    all_ids.append(sample_id)
                    all_scales.append(scale_factor)
                    all_depths.append(normalized_depths)
        
        except Exception as e:
            log_message(f"Error reading batch {batch}: {e}")
            continue
        
        log_message(f"Read batch {batch}")
    
    if not all_depths:
        log_message("Error: No samples successfully loaded!")
        return [], [], np.array([])
    
    depths_matrix = np.array(all_depths)
    log_message(f"Read {len(all_ids)} individuals; normalizing by region")
    
    return all_ids, all_scales, depths_matrix

def normalize_by_region(depths_matrix):
    """Normalize each region across individuals"""
    n_samples, n_regions = depths_matrix.shape
    
    normalized_depths = np.zeros_like(depths_matrix)
    region_stats = []
    
    for r in range(n_regions):
        region_depths = depths_matrix[:, r]
        
        mu = np.mean(region_depths)
        variance = np.var(region_depths, ddof=1) if n_samples > 1 else 0
        
        region_stats.append({
            'mean': mu,
            'variance': variance,
            'sigma2_ratio': 100 * variance / mu if mu > 0 else 0
        })
        
        # Normalize: (x - mu) / sqrt(mu) - matching C++ logic
        if mu > 0:
            inv_root_mean = 1.0 / math.sqrt(mu)
            normalized_depths[:, r] = (region_depths - mu) * inv_root_mean
        else:
            normalized_depths[:, r] = region_depths - mu
    
    log_message("Normalized by region")
    return normalized_depths, region_stats

def select_high_variance_regions(region_stats, percentile=0.9):
    """Select regions with high variance for final output"""
    ratio_mult = 100
    sigma2_ratios = [stats['sigma2_ratio'] for stats in region_stats]
    
    if not sigma2_ratios:
        return [], 1.0
    
    sigma2_ratios_sorted = sorted(sigma2_ratios)
    
    # Select regions above 90th percentile
    threshold_idx = int(percentile * len(sigma2_ratios_sorted))
    if threshold_idx >= len(sigma2_ratios_sorted):
        threshold_idx = len(sigma2_ratios_sorted) - 1
    sigma2_ratio_min = sigma2_ratios_sorted[threshold_idx]
    
    want_regions = [ratio > sigma2_ratio_min for ratio in sigma2_ratios]
    num_want = sum(want_regions)
    
    log_message(f"Restricting to {num_want} regions with sigma2ratio > {sigma2_ratio_min:.3f}")
    
    # Calculate scaling factor based on median sigma2ratio
    median_idx = len(sigma2_ratios_sorted) // 2
    sigma2_ratio_median = sigma2_ratios_sorted[median_idx] if sigma2_ratios_sorted else 1.0
    scale_factor = 1.0 / math.sqrt(sigma2_ratio_median / ratio_mult) if sigma2_ratio_median > 0 else 1.0
    
    log_message(f"Rescaling to approximate z-scores based on median sigma2ratio = {sigma2_ratio_median:.3f}")
    
    return want_regions, scale_factor

def write_output(output_file, sample_ids, scales, normalized_depths, region_stats, want_regions, scale_factor):
    """Write normalized output to gzipped file"""
    log_message(f"Writing output to {output_file}")
    
    n_samples = len(sample_ids)
    selected_regions = [i for i, want in enumerate(want_regions) if want]
    n_regions = len(selected_regions)
    
    with gzip.open(output_file, 'wt') as f:
        # Write region means
        f.write(f"{n_samples}\t{n_regions}")
        for r in selected_regions:
            f.write(f"\t{region_stats[r]['mean']:.3f}")
        f.write("\n")
        
        # Write sigma2 ratios
        f.write(f"{n_samples}\t{n_regions}")
        for r in selected_regions:
            f.write(f"\t{region_stats[r]['sigma2_ratio']:.3f}")
        f.write("\n")
        
        # Write sample data
        for i, sample_id in enumerate(sample_ids):
            f.write(f"{sample_id}\t{scales[i]*0.01:.2f}")
            for r in selected_regions:
                scaled_depth = scale_factor * normalized_depths[i, r]
                f.write(f"\t{scaled_depth:.2f}")
            f.write("\n")

def main():
    if len(sys.argv) != 6:
        print("ERROR: 5 arguments required")
        print("- arg1: prefix of mosdepth input (no more than 170 characters), used in <prefix>_batch_<batchnumber>.txt.gz")
        print("- arg2: bed file path e.g. /path/to/repeat_mask_list.hg38.ucsc_bed")
        print("- arg3: example input e.g. /path/to/name_regions.bed.gz")
        print("- arg4: N_sample(int)")
        print("- arg5: output path e.g. /path/to/ID_scale_zdepths.txt.gz")
        return 1
    
    mosdepth_prefix = sys.argv[1]
    bed_source = sys.argv[2]
    example_output = sys.argv[3]
    max_n_samples = int(sys.argv[4])
    output_path = sys.argv[5]
    
    start_time = time.time()
    
    # Step 1: Determine which regions to extract based on first 10 batches
    extract_regions, mean_depths, n_analyzed = determine_regions_to_extract(mosdepth_prefix)
    
    if not extract_regions or not any(extract_regions):
        log_message("Error: No regions meet extraction criteria!")
        return 1
    
    # Step 2: Load repeat mask and exclude overlapping regions
    overlaps_vntr = load_repeat_mask(bed_source)
    extract_regions = update_extraction_list(example_output, extract_regions, overlaps_vntr)
    
    if not any(extract_regions):
        log_message("Error: No regions remain after VNTR filtering!")
        return 1
    
    # Step 3: Load all sample data and normalize by individual
    sample_ids, scales, depths_matrix = load_all_samples(mosdepth_prefix, extract_regions, max_n_samples)
    
    if len(sample_ids) == 0:
        log_message("Error: No samples successfully loaded!")
        return 1
    
    # Step 4: Normalize by region across individuals
    normalized_depths, region_stats = normalize_by_region(depths_matrix)
    
    # Step 5: Select high-variance regions
    want_regions, scale_factor = select_high_variance_regions(region_stats)
    
    # Step 6: Write output
    write_output(output_path, sample_ids, scales, normalized_depths, region_stats, want_regions, scale_factor)
    
    total_time = time.time() - start_time
    log_message(f"Complete! Total time: {total_time:.1f} seconds")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())