#!/usr/bin/env python3
"""
Compare mosdepth coverage for two LPA regions (hg38 vs samtools coords).
Appends results (mean coverage values, difference, percent difference) to a table.
"""

import sys
import gzip
import time
from pathlib import Path
import argparse


def setup_logging(log_file):
    """Redirect stdout to a log file"""
    sys.stdout = open(log_file, 'a')
    log_message(f"Logging to {log_file}")


def log_message(msg):
    """Print timestamped log message"""
    print(f"[{time.strftime('%H:%M:%S')}] {msg}")


def read_gzipped_file(filepath):
    """Open gzipped or plain file"""
    if str(filepath).endswith('.gz'):
        return gzip.open(filepath, 'rt')
    else:
        return open(filepath, 'r')


def extract_depths(filepath, chromosome, start, end):
    """Extract depth values overlapping a region from mosdepth .regions.bed.gz"""
    depths = []
    with read_gzipped_file(filepath) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split('\t')
            if len(parts) < 4:
                continue

            chr_str, line_start, line_end, depth = parts[0], int(parts[1]), int(parts[2]), float(parts[3])
            if chr_str == chromosome and line_start <= end and line_end >= start:
                depths.append(depth)

    return depths


def parse_args():
    parser = argparse.ArgumentParser(description="Compare mosdepth coverage in two regions.")
    parser.add_argument("mosdepth_bed", type=str, help="Path to mosdepth .regions.bed.gz file")
    parser.add_argument("chr", type=str, help="Chromosome")
    parser.add_argument("hg38_start", type=int, help="HG38 start")
    parser.add_argument("hg38_end", type=int, help="HG38 end")
    parser.add_argument("sam_start", type=int, help="Samtools start")
    parser.add_argument("sam_end", type=int, help="Samtools end")
    parser.add_argument("output_file", type=str, help="Output table file")
    parser.add_argument("log_file", type=str, help="Log file")
    return parser.parse_args()


def main():
    args = parse_args()
    setup_logging(args.log_file)

    log_message(f"Analyzing {args.mosdepth_bed}")
    hg38_depths = extract_depths(args.mosdepth_bed, f"chr{args.chr}", args.hg38_start, args.hg38_end)
    sam_depths = extract_depths(args.mosdepth_bed, f"chr{args.chr}", args.sam_start, args.sam_end)

    mean1 = sum(hg38_depths) / len(hg38_depths) if hg38_depths else 0
    mean2 = sum(sam_depths) / len(sam_depths) if sam_depths else 0
    diff = mean1 - mean2
    pct_diff = (diff / mean1 * 100) if mean1 != 0 else 0

    log_message(f"HG38 mean depth: {mean1:.3f} (n={len(hg38_depths)})")
    log_message(f"Samtools mean depth: {mean2:.3f} (n={len(sam_depths)})")
    log_message(f"Difference: {diff:.3f}, Percent diff: {pct_diff:.2f}%")

    with open(args.output_file, 'a') as out_f:
        out_f.write(
            f"{Path(args.mosdepth_bed).stem}\t{mean1:.3f}\t{mean2:.3f}\t{diff:.3f}\t{pct_diff:.2f}\n"
        )


if __name__ == "__main__":
    main()
