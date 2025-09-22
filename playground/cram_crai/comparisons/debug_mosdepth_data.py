#!/usr/bin/env python3
"""
Debug script to analyze mosdepth data specifically for the LPA gene using logging
"""

import sys
import gzip
import glob
from pathlib import Path
import time

def setup_logging(log_file):
    """Setup logging to file and console"""
    sys.stdout = open(log_file, 'a')
    log_message(f"Logging to {log_file}")

def log_message(msg):
    """Print timestamped log message"""
    print(f"[{time.strftime('%H:%M:%S')}] {msg}")

def read_gzipped_file(filepath):
    """Read gzipped or regular file"""
    if str(filepath).endswith('.gz'):
        return gzip.open(filepath, 'rt')
    else:
        return open(filepath, 'r')

def analyze_sample_file(
    filepath,
    chromosome='chr6',
    start=160605062,
    end=160647661
):
    """Analyze a single regions file for a gene only"""
    log_message(f"\nAnalyzing: {filepath}")
    
    lpa_regions = 0
    lpa_depths = []
    
    try:
        with read_gzipped_file(filepath) as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                
                parts = line.strip().split('\t')
                if len(parts) < 4:
                    continue
                
                chr_str = parts[0]
                line_start = int(parts[1])
                line_end = int(parts[2])
                depth = float(parts[3])
                
                # Bin overlaps gene if line_start <= end and line_end >= start
                if chr_str == chromosome and line_start <= end and line_end >= start:
                    lpa_regions += 1
                    lpa_depths.append(depth)
                    
                    # Show first few lines as examples
                    if len(lpa_depths) <= 5:
                        log_message(f"  LPA Line: {chr_str}:{line_start}-{line_end} depth={depth}")
    
    except Exception as e:
        log_message(f"Error reading {filepath}: {e}")
        return
    
    log_message(f"{chromosome} gene regions: {lpa_regions}")
    if lpa_depths:
        log_message(f"{chromosome} depths - Min: {min(lpa_depths):.3f}, Max: {max(lpa_depths):.3f}, Mean: {sum(lpa_depths)/len(lpa_depths):.3f}")
        log_message(f"Sample {chromosome} depths: {lpa_depths[:10]}")
    else:
        log_message("No gene regions found!")


def main():
    if len(sys.argv) <= 3:
        print("Usage: python debug_mosdepth_data.py <mosdepth_dir> <log_file>")
        return 1
    
    mosdepth_dir = sys.argv[1]
    log_file = sys.argv[2]
    chromosome=sys.argv[3] if len(sys.argv) > 3 else 'chr6'
    start = int(sys.argv[4]) if len(sys.argv) > 4 else 160605062
    end   = int(sys.argv[5]) if len(sys.argv) > 5 else 160647661


    setup_logging(log_file)
    
    # Find regions files
    pattern = str(Path(mosdepth_dir) / "*.regions.bed.gz")
    files = glob.glob(pattern)
    
    if not files:
        log_message(f"No .regions.bed.gz files found in {mosdepth_dir}")
        return 1
    
    log_message(f"Found {len(files)} files")
    
    # Analyze first few files
    for i, file_path in enumerate(sorted(files)[:3]):
        analyze_sample_file(file_path, chromosome, start, end)
        if i >= 2:  # Stop after 3 files
            break

if __name__ == "__main__":
    sys.exit(main())
