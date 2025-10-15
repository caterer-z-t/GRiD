#!/usr/bin/env python3
"""
Compute Diploid Copy Number for LPA KIV2 Repeats

This script normalizes read counts from realignment using neighbor-based normalization.
For each exon type, it finds the top N neighbors and normalizes counts using both
individual scale and neighbor mean scale.
"""

import argparse
import gzip
import sys
from pathlib import Path


def load_count_results(count_file):
    """
    Load read counts from realignment output.
    
    Format: sample_id\t1B_KIV3\t1B_KIV2\t1B_tied\t1A
    
    Returns dict mapping sample_id to dict of counts
    """
    counts = {}
    
    with open(count_file, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) != 5:
                continue
            
            sample_id = fields[0]
            try:
                counts[sample_id] = {
                    '1B_KIV3': int(fields[1]),
                    '1B_KIV2': int(fields[2]),
                    '1B_tied': int(fields[3]),
                    '1A': int(fields[4])
                }
            except ValueError:
                print(f"WARNING: Skipping malformed line: {line.strip()}", file=sys.stderr)
                continue
    
    return counts


def load_neighbor_results(neighbor_file):
    """
    Load neighbor information from find_neighbors output.
    
    Format: sample_id\tscale\tneighbor1\tneighbor1_scale\tneighbor1_distance\t...
    
    Returns dict mapping sample_id to (scale, [(neighbor_id, neighbor_scale, distance), ...])
    """
    neighbors = {}
    
    # Handle both gzipped and plain text files
    open_func = gzip.open if neighbor_file.endswith('.gz') else open
    mode = 'rt' if neighbor_file.endswith('.gz') else 'r'
    
    with open_func(neighbor_file, mode) as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) < 2:
                continue
            
            sample_id = fields[0]
            try:
                scale = float(fields[1])
            except ValueError:
                continue
            
            # Parse neighbors (triplets: neighbor_id, scale, distance)
            neighbor_list = []
            for i in range(2, len(fields), 3):
                if i + 2 < len(fields):
                    try:
                        neighbor_id = fields[i]
                        neighbor_scale = float(fields[i + 1])
                        distance = float(fields[i + 2])
                        neighbor_list.append((neighbor_id, neighbor_scale, distance))
                    except ValueError:
                        continue
            
            neighbors[sample_id] = (scale, neighbor_list)
    
    return neighbors


def get_exon_count(counts_dict, exon_type):
    """
    Get the appropriate count for a given exon type.
    
    Exon types:
    - 1B_KIV3: just the 1B_KIV3 count
    - 1B_notKIV3: 1B_KIV2 + 1B_tied
    - 1B: 1B_KIV3 + 1B_KIV2 + 1B_tied
    - 1A: just the 1A count
    """
    if exon_type == '1B_KIV3':
        return counts_dict.get('1B_KIV3', 0)
    elif exon_type == '1B_notKIV3':
        return counts_dict.get('1B_KIV2', 0) + counts_dict.get('1B_tied', 0)
    elif exon_type == '1B':
        return (counts_dict.get('1B_KIV3', 0) + 
                counts_dict.get('1B_KIV2', 0) + 
                counts_dict.get('1B_tied', 0))
    elif exon_type == '1A':
        return counts_dict.get('1A', 0)
    else:
        raise ValueError(f"Unknown exon type: {exon_type}")


def compute_diploid_cn(counts, neighbors, exon_type, n_neighbors=200):
    """
    Compute diploid copy number for all samples for a given exon type.
    
    Formula: dipCN = (sample_count / sample_scale) / (mean_neighbor_normalized_count)
    where mean_neighbor_normalized_count = mean(neighbor_count / neighbor_scale)
    
    Returns dict mapping sample_id to diploid copy number
    """
    results = {}
    
    for sample_id, (sample_scale, neighbor_list) in neighbors.items():
        # Get count for this sample
        if sample_id not in counts:
            continue
        
        sample_count = get_exon_count(counts[sample_id], exon_type)
        
        if sample_count == 0:
            # Skip samples with zero counts
            continue
        
        # Compute mean normalized count across top N neighbors
        neighbor_sum = 0.0
        neighbor_num = 0
        
        for i, (neighbor_id, neighbor_scale, distance) in enumerate(neighbor_list[:n_neighbors]):
            if neighbor_id not in counts:
                continue
            
            neighbor_count = get_exon_count(counts[neighbor_id], exon_type)
            
            if neighbor_count > 0 and neighbor_scale > 0:
                neighbor_sum += neighbor_count / neighbor_scale
                neighbor_num += 1
        
        # Calculate diploid copy number
        if neighbor_num > 0 and sample_scale > 0:
            mean_neighbor_normalized = neighbor_sum / neighbor_num
            if mean_neighbor_normalized > 0:
                dip_cn = (sample_count / sample_scale) / mean_neighbor_normalized
                results[sample_id] = dip_cn
    
    return results


def write_output(results, output_file):
    """Write diploid copy numbers to output file."""
    with open(output_file, 'w') as f:
        f.write("ID\tdipCN\n")
        for sample_id, dip_cn in sorted(results.items()):
            f.write(f"{sample_id}\t{dip_cn:.6f}\n")


def main():
    parser = argparse.ArgumentParser(
        description='Compute diploid copy numbers for LPA KIV2 repeats'
    )
    parser.add_argument(
        '--count_file',
        required=True,
        help='Realignment output file with read counts'
    )
    parser.add_argument(
        '--neighbor_file',
        required=True,
        help='Find neighbors output file (can be gzipped)'
    )
    parser.add_argument(
        '--output_prefix',
        required=True,
        help='Output file prefix (will create 4 files: .exon1A.dipCN.txt, etc.)'
    )
    parser.add_argument(
        '--n_neighbors',
        type=int,
        default=200,
        help='Number of top neighbors to use (default: 200)'
    )
    
    args = parser.parse_args()
    
    # Load input files
    print(f"Loading count results from {args.count_file}...")
    counts = load_count_results(args.count_file)
    print(f"  Loaded counts for {len(counts)} samples")
    
    print(f"Loading neighbor results from {args.neighbor_file}...")
    neighbors = load_neighbor_results(args.neighbor_file)
    print(f"  Loaded neighbor info for {len(neighbors)} samples")
    
    # Process each exon type
    exon_types = ['1B_KIV3', '1B_notKIV3', '1B', '1A']
    
    for exon_type in exon_types:
        print(f"\nProcessing exon type: {exon_type}")
        
        # Compute diploid copy numbers
        results = compute_diploid_cn(counts, neighbors, exon_type, args.n_neighbors)
        
        # Write output
        output_file = f"{args.output_prefix}.exon{exon_type}.dipCN.txt"
        write_output(results, output_file)
        
        print(f"  Computed dipCN for {len(results)} samples")
        print(f"  Output written to: {output_file}")
    
    print("\nAll exon types processed successfully!")


if __name__ == '__main__':
    main()