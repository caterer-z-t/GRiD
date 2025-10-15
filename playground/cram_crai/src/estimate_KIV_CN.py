#!/usr/bin/env python3
"""
Estimate LPA KIV2 Copy Number

This script computes the final KIV2 copy number estimates using the formula:
    diploid_estimate = 34.9 × exon1A + 5.2 × exon1B - 1
    haploid_estimate = diploid_estimate / 2

The coefficients were empirically determined to best predict copy number
from the normalized read depth values.
"""

import argparse
import pandas as pd
import sys
from pathlib import Path


def load_dipcn_file(file_path, exon_name):
    """
    Load diploid copy number file.
    
    Expected format:
        ID    dipCN
        sample1    15.234
        sample2    18.456
    
    Returns DataFrame with index=ID and column=exon_name
    """
    try:
        # Try tab-separated first (our output format)
        df = pd.read_csv(file_path, sep='\t', index_col=0)
    except:
        try:
            # Try space-separated (original format)
            df = pd.read_csv(file_path, sep=' ', index_col=0)
        except Exception as e:
            print(f"ERROR: Failed to load {file_path}: {e}", file=sys.stderr)
            sys.exit(1)
    
    # Rename column to exon name
    df.rename(columns={'dipCN': exon_name}, inplace=True)
    
    return df


def compute_kiv_estimates(exon1a_file, exon1b_file):
    """
    Compute KIV2 copy number estimates from exon1A and exon1B diploid copy numbers.
    
    Formula:
        diploid_estimate = 34.9 × exon1A + 5.2 × exon1B - 1
        haploid_estimate = diploid_estimate / 2
    
    Returns DataFrame with columns: exon1A, exon1B, dip_estimate, estimate
    """
    # Load diploid copy numbers
    print(f"Loading exon1A diploid copy numbers from {exon1a_file}")
    exon1a = load_dipcn_file(exon1a_file, 'exon1A')
    print(f"  Loaded {len(exon1a)} samples")
    
    print(f"Loading exon1B diploid copy numbers from {exon1b_file}")
    exon1b = load_dipcn_file(exon1b_file, 'exon1B')
    print(f"  Loaded {len(exon1b)} samples")
    
    # Merge on sample IDs (inner join - only samples in both files)
    print("\nMerging exon1A and exon1B data...")
    combined = pd.concat([exon1a, exon1b], join='inner', axis=1)
    print(f"  {len(combined)} samples have both exon1A and exon1B data")
    
    if len(combined) == 0:
        print("ERROR: No overlapping samples between exon1A and exon1B files!", file=sys.stderr)
        sys.exit(1)
    
    # Compute estimates using the empirical formula
    print("\nComputing KIV2 copy number estimates...")
    combined['dip_estimate'] = 34.9 * combined['exon1A'] + 5.2 * combined['exon1B'] - 1
    combined['estimate'] = combined['dip_estimate'] / 2
    
    # Summary statistics
    print("\nSummary statistics:")
    print(f"  Diploid estimates - Mean: {combined['dip_estimate'].mean():.2f}, "
          f"Median: {combined['dip_estimate'].median():.2f}, "
          f"Std: {combined['dip_estimate'].std():.2f}")
    print(f"  Haploid estimates - Mean: {combined['estimate'].mean():.2f}, "
          f"Median: {combined['estimate'].median():.2f}, "
          f"Std: {combined['estimate'].std():.2f}")
    print(f"  Range: {combined['estimate'].min():.2f} - {combined['estimate'].max():.2f}")
    
    return combined


def main():
    parser = argparse.ArgumentParser(
        description='Estimate LPA KIV2 copy numbers from diploid copy number values',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Formula:
    diploid_estimate = 34.9 × exon1A + 5.2 × exon1B - 1
    haploid_estimate = diploid_estimate / 2

Output columns:
    exon1A          - Normalized diploid copy number for exon1A
    exon1B          - Normalized diploid copy number for exon1B
    dip_estimate    - Diploid KIV2 copy number estimate
    estimate        - Haploid KIV2 copy number estimate
        """
    )
    
    parser.add_argument(
        '--exon1a',
        required=True,
        help='Path to exon1A diploid CN file (*.exon1A.dipCN.txt)'
    )
    parser.add_argument(
        '--exon1b',
        required=True,
        help='Path to exon1B diploid CN file (*.exon1B.dipCN.txt)'
    )
    parser.add_argument(
        '--output',
        required=True,
        help='Output file path for KIV2 copy number estimates'
    )
    parser.add_argument(
        '--format',
        choices=['csv', 'tsv', 'txt'],
        default='tsv',
        help='Output format (default: tsv)'
    )
    
    args = parser.parse_args()
    
    # Validate input files exist
    if not Path(args.exon1a).exists():
        print(f"ERROR: exon1A file not found: {args.exon1a}", file=sys.stderr)
        sys.exit(1)
    
    if not Path(args.exon1b).exists():
        print(f"ERROR: exon1B file not found: {args.exon1b}", file=sys.stderr)
        sys.exit(1)
    
    # Compute estimates
    results = compute_kiv_estimates(args.exon1a, args.exon1b)
    
    # Create output directory if needed
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Write output
    print(f"\nWriting results to {args.output}")
    if args.format == 'csv':
        results.to_csv(args.output, index=True)
    elif args.format == 'tsv' or args.format == 'txt':
        results.to_csv(args.output, index=True, sep='\t')
    
    print(f"Successfully wrote {len(results)} samples to output file")
    
    # Display first few rows
    print("\nFirst 10 samples:")
    print(results.head(10).to_string())


if __name__ == '__main__':
    main()