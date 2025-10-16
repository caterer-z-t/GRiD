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


def normalize_sample_id(sample_id):
    """
    Normalize sample IDs by removing common file extensions and suffixes.
    
    Examples:
        NWD278973.b38.irc.v1_subset → NWD278973
        NWD278973.cram → NWD278973
        NWD278973 → NWD278973
    """
    # Strip whitespace first
    sample_id = sample_id.strip()
    
    # Remove .b38.irc.v1_subset pattern (before other extensions)
    if '.b38.irc.v1_subset' in sample_id:
        sample_id = sample_id.replace('.b38.irc.v1_subset', '')
    
    # Remove common CRAM/BAM file patterns
    if sample_id.endswith('.cram'):
        sample_id = sample_id[:-5]
    elif sample_id.endswith('.bam'):
        sample_id = sample_id[:-4]
    
    return sample_id.strip()


def load_count_results(count_file):
    """
    Load read counts from realignment output.
    
    Format: sample_id\t1B_KIV3\t1B_KIV2\t1B_tied\t1A
    
    Returns dict mapping sample_id to dict of counts
    """
    counts = {}
    
    with open(count_file, 'r') as f:
        for i, line in enumerate(f):
            fields = line.strip().split('\t')
            if len(fields) != 5:
                continue
            
            original_id = fields[0]
            sample_id = normalize_sample_id(original_id)
            
            # Debug first few
            if i < 3:
                print(f"  DEBUG count: '{original_id}' -> '{sample_id}'", file=sys.stderr)
            
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
        for i, line in enumerate(f):
            fields = line.strip().split('\t')
            if len(fields) < 2:
                continue
            
            original_id = fields[0]
            sample_id = normalize_sample_id(original_id)
            
            # Debug first few
            if i < 3:
                print(f"  DEBUG neighbor: '{original_id}' -> '{sample_id}'", file=sys.stderr)
            
            try:
                scale = float(fields[1])
            except ValueError:
                continue
            
            # Parse neighbors (triplets: neighbor_id, scale, distance)
            neighbor_list = []
            for j in range(2, len(fields), 3):
                if j + 2 < len(fields):
                    try:
                        neighbor_id = normalize_sample_id(fields[j])
                        neighbor_scale = float(fields[j + 1])
                        distance = float(fields[j + 2])
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
    
    # Debug: Check for sample ID overlap
    count_ids = set(counts.keys())
    neighbor_ids = set(neighbors.keys())
    overlap = count_ids & neighbor_ids
    
    print(f"\nSample ID debugging:")
    print(f"  Samples in count file: {len(count_ids)}")
    print(f"  Samples in neighbor file: {len(neighbor_ids)}")
    print(f"  Overlapping samples: {len(overlap)}")
    
    if len(overlap) == 0:
        print("\nERROR: No overlapping sample IDs found!")
        print("\nFirst 5 sample IDs from count file:")
        for i, sample_id in enumerate(list(count_ids)[:5]):
            print(f"    '{sample_id}'")
        print("\nFirst 5 sample IDs from neighbor file:")
        for i, sample_id in enumerate(list(neighbor_ids)[:5]):
            print(f"    '{sample_id}'")
        print("\nSample IDs don't match! Check if one has file extensions or different naming.")
        sys.exit(1)
    elif len(overlap) < len(count_ids) * 0.5:
        print(f"\nWARNING: Only {len(overlap)/len(count_ids)*100:.1f}% of samples overlap!")
        print("  Some samples in count file may be missing from neighbor file")
    
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