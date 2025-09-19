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
import pandas as pd
from pathlib import Path
import time
import glob
from collections import defaultdict

def log_message(msg):
    """Print timestamped log message"""
    print(f"[{time.strftime('%H:%M:%S')}] {msg}")

def read_gzipped_file(filepath):
    """Read gzipped or regular file"""
    if str(filepath).endswith('.gz'):
        return gzip.open(filepath, 'rt')
    else:
        return open(filepath, 'r')

def parse_chromosome(chr_str):
    """Parse chromosome string to integer"""
    if chr_str.startswith('chr'):
        chr_str = chr_str[3:]
    try:
        return int(chr_str)
    except ValueError:
        return None  # Skip non-numeric chromosomes

def find_regions_files(mosdepth_dir):
    """Find all .regions.bed.gz files and extract sample IDs"""
    regions_files = []
    sample_ids = []
    
    # Look for .regions.bed.gz files
    pattern = str(Path(mosdepth_dir) / "*.regions.bed.gz")
    files = glob.glob(pattern)
    
    for file_path in sorted(files):
        filename = Path(file_path).name
        # Extract sample ID from filename (everything before the first dot or specific suffix)
        if filename.endswith('.regions.bed.gz'):
            # For your format: NWD997185.b38.irc.v1_subset_LPA.regions.bed.gz
            sample_id = filename.split('.')[0]  # Gets NWD997185
            regions_files.append(file_path)
            sample_ids.append(sample_id)
    
    log_message(f"Found {len(regions_files)} regions files")
    return regions_files, sample_ids

def read_individual_regions_file(regions_file, target_chr=6):
    """Read a single individual's .regions.bed.gz file and return depths for target chromosome"""
    depths = []
    regions_info = []
    
    try:
        with read_gzipped_file(regions_file) as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                
                parts = line.strip().split('\t')
                if len(parts) < 4:
                    continue
                
                chr_num = parse_chromosome(parts[0])
                if chr_num != target_chr:
                    continue
                
                start = int(parts[1])
                end = int(parts[2])
                depth = float(parts[3])
                
                depths.append(depth)
                regions_info.append((chr_num, start, end))
        
        return depths, regions_info
    except FileNotFoundError:
        log_message(f"Warning: Regions file {regions_file} not found, skipping")
        return [], []

def determine_regions_to_extract(regions_files, sample_ids, max_samples_for_analysis=250):
    """
    Analyze subset of samples to determine which regions to extract
    based on mean depth criteria (20-100x coverage)
    """
    log_message("Analyzing samples to determine regions to extract...")
    
    # Use first N samples for analysis (or all if fewer than max)
    analysis_files = regions_files[:min(len(regions_files), max_samples_for_analysis)]
    analysis_ids = sample_ids[:len(analysis_files)]
    
    all_depths = []
    region_coords = []
    total_samples = 0
    
    for i, (regions_file, sample_id) in enumerate(zip(analysis_files, analysis_ids)):
        depths, regions_info = read_individual_regions_file(regions_file)
        
        if not depths:
            continue
        
        if total_samples == 0:
            # Initialize with first sample
            all_depths = [0.0] * len(depths)
            region_coords = regions_info
        
        if len(depths) != len(all_depths):
            log_message(f"Warning: Sample {sample_id} has different number of regions ({len(depths)} vs {len(all_depths)})")
            continue
        
        # Add depths (don't convert yet - let's see raw values first)
        for r, depth in enumerate(depths):
            all_depths[r] += depth
        
        total_samples += 1
        
        if (i + 1) % 50 == 0:
            log_message(f"Processed {i + 1} samples for region analysis")
    
    # Calculate mean depths and determine extraction criteria
    if total_samples > 0:
        mean_depths = [d / total_samples for d in all_depths]
        
        # Debug: show some statistics about the raw depths
        raw_depths_sample = mean_depths[:20]  # First 20 regions
        log_message(f"Sample of raw mean depths: {[f'{d:.2f}' for d in raw_depths_sample]}")
        log_message(f"Min/Max/Median raw depth: {min(mean_depths):.2f}/{max(mean_depths):.2f}/{sorted(mean_depths)[len(mean_depths)//2]:.2f}")
        
        # Try different thresholds to see what might work
        for min_thresh, max_thresh in [(0.2, 10), (2, 100), (20, 1000), (0.02, 1)]:
            candidate_regions = [min_thresh <= d <= max_thresh for d in mean_depths]
            log_message(f"Regions with depth {min_thresh}-{max_thresh}x: {sum(candidate_regions)}")
        
        # Use the original threshold but remove the 0.01 conversion
        extract_regions = [20 <= mean_depth <= 100 for mean_depth in mean_depths]
        
    else:
        mean_depths = []
        extract_regions = []
    
    num_extract = sum(extract_regions)
    log_message(f"Analyzed {total_samples} individuals")
    log_message(f"Extracting {num_extract} / {len(all_depths)} regions with depth 20-100x")
    
    return extract_regions, mean_depths, region_coords

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
        
        log_message(f"Read {num_vntrs} VNTRs spanning {total_length/1e6:.2f} Mb")
        return overlaps_vntr
    
    except FileNotFoundError:
        log_message(f"Warning: Repeat mask file {bed_file} not found")
        return defaultdict(lambda: defaultdict(bool))

def update_extraction_list(region_coords, extract_regions, overlaps_vntr):
    """Update extraction list to exclude regions overlapping VNTRs"""
    log_message("Updating extraction list to exclude VNTR overlaps...")
    
    regions_overlap = 0
    extract_overlap = 0
    
    for r, (chr_num, start, end) in enumerate(region_coords):
        if r >= len(extract_regions):
            break
        
        # Check if region overlaps VNTR
        if overlaps_vntr[chr_num][start // 1000]:
            regions_overlap += 1
            if extract_regions[r]:
                extract_overlap += 1
                extract_regions[r] = False
    
    num_extract = sum(extract_regions)
    log_message(f"Excluding {regions_overlap} / {len(region_coords)} regions overlapping VNTRs")
    log_message(f"Excluded {extract_overlap} in extract set; {num_extract} left")
    
    return extract_regions

def load_all_samples(regions_files, sample_ids, extract_regions, max_samples):
    """Load depth data for all samples"""
    log_message("Loading all sample data...")
    
    num_extract = sum(extract_regions)
    if num_extract == 0:
        log_message("Error: No regions selected for extraction!")
        return [], [], np.array([])
    
    # Limit to max_samples
    files_to_process = regions_files[:max_samples]
    ids_to_process = sample_ids[:max_samples]
    
    all_ids = []
    all_scales = []
    all_depths = np.zeros((len(files_to_process), num_extract), dtype=np.float32)
    
    successful_samples = 0
    
    for i, (regions_file, sample_id) in enumerate(zip(files_to_process, ids_to_process)):
        depths, _ = read_individual_regions_file(regions_file)
        
        if not depths or len(depths) != len(extract_regions):
            log_message(f"Warning: Skipping {sample_id} - invalid depth data")
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
        all_depths[successful_samples, :] = normalized_depths
        
        successful_samples += 1
        
        if (successful_samples) % 100 == 0:
            log_message(f"Processed {successful_samples} samples")
    
    log_message(f"Successfully loaded {successful_samples} individuals; normalizing by region")
    
    return all_ids, all_scales, all_depths[:successful_samples, :]

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
        
        # Normalize: (x - mu) / sqrt(mu)
        if mu > 0:
            inv_root_mean = 1.0 / math.sqrt(mu)
            normalized_depths[:, r] = (region_depths - mu) * inv_root_mean
        else:
            normalized_depths[:, r] = region_depths - mu
    
    return normalized_depths, region_stats

def select_high_variance_regions(region_stats, percentile=0.9):
    """Select regions with high variance for final output"""
    ratio_mult = 100
    sigma2_ratios = [stats['sigma2_ratio'] for stats in region_stats]
    
    if not sigma2_ratios:
        return [], 1.0
    
    sigma2_ratios_sorted = sorted(sigma2_ratios)
    
    threshold_idx = int(percentile * len(sigma2_ratios_sorted))
    sigma2_ratio_min = sigma2_ratios_sorted[threshold_idx] if threshold_idx < len(sigma2_ratios_sorted) else 0
    
    want_regions = [ratio >= sigma2_ratio_min for ratio in sigma2_ratios]
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
    if len(sys.argv) != 5:
        print("ERROR: 4 arguments required")
        print("- arg1: mosdepth directory containing .regions.bed.gz files")
        print("- arg2: bed file path e.g. /path/to/repeat_mask_list.bed")
        print("- arg3: N_sample(int) - maximum number of samples to process")
        print("- arg4: output path e.g. /path/to/ID_scale_zdepths.txt.gz")
        return 1
    
    mosdepth_dir = sys.argv[1]
    bed_source = sys.argv[2]
    max_samples = int(sys.argv[3])
    output_path = sys.argv[4]
    
    start_time = time.time()
    
    # Step 1: Find all regions files
    regions_files, sample_ids = find_regions_files(mosdepth_dir)
    
    if not regions_files:
        log_message("Error: No .regions.bed.gz files found!")
        return 1
    
    # Step 2: Determine which regions to extract
    extract_regions, mean_depths, region_coords = determine_regions_to_extract(regions_files, sample_ids)
    
    if not extract_regions or not any(extract_regions):
        log_message("Error: No regions meet extraction criteria!")
        return 1
    
    # Step 3: Load repeat mask and exclude overlapping regions
    overlaps_vntr = load_repeat_mask(bed_source)
    extract_regions = update_extraction_list(region_coords, extract_regions, overlaps_vntr)
    
    if not any(extract_regions):
        log_message("Error: No regions remain after VNTR filtering!")
        return 1
    
    # Step 4: Load all sample data and normalize by individual
    sample_ids_final, scales, depths_matrix = load_all_samples(regions_files, sample_ids, extract_regions, max_samples)
    
    if len(sample_ids_final) == 0:
        log_message("Error: No samples successfully loaded!")
        return 1
    
    # Step 5: Normalize by region across individuals
    normalized_depths, region_stats = normalize_by_region(depths_matrix)
    
    # Step 6: Select high-variance regions
    want_regions, scale_factor = select_high_variance_regions(region_stats)
    
    # Step 7: Write output
    write_output(output_path, sample_ids_final, scales, normalized_depths, region_stats, want_regions, scale_factor)
    
    total_time = time.time() - start_time
    log_message(f"Complete! Total time: {total_time:.1f} seconds")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())