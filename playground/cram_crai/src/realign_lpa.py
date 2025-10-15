#!/usr/bin/env python3
"""
LPA Realignment Script - Python version
Converts C++ realignment logic to Python for easier maintenance and integration.

This script processes samtools output to classify LPA variants based on read depth
at specific positions in the LPA gene region.
"""

import sys
import argparse
from collections import defaultdict
import re

# Constants
MIN_QUAL = ord(':')  # ASCII value for ':'
BAD_SCORE = 9999

# Positions of the 7 KIV2 repeats in hg19/hg37
STARTS = [161032593, 161038148, 161043694, 161049238, 161054784, 161060331, 161065878]
REF_OFFSET = 161032032  # Reference start position


def load_reference(ref_fasta_path):
    """Load reference sequence from FASTA file."""
    with open(ref_fasta_path, 'r') as f:
        lines = f.readlines()
        # Skip header line(s) and concatenate sequence
        ref_seq = ''.join(line.strip() for line in lines if not line.startswith('>'))
    return ref_seq


def parse_cigar(cigar):
    """Parse CIGAR string into list of (length, operation) tuples."""
    pattern = re.compile(r'(\d+)([MIDNSHP=X])')
    return [(int(length), op) for length, op in pattern.findall(cigar)]


def compute_scores(pos, cigar, seq, qual, ref):
    """
    Compute mismatch scores for each of the 7 repeat positions.
    
    Returns array of 7 scores, one for each potential repeat alignment.
    Lower score = better match to that repeat position.
    """
    scores = [0] * 7
    
    # Find which repeat this read likely aligns to based on position
    rep_bwa = -1
    for i, start in enumerate(STARTS):
        if start - 500 < pos < start + 500:
            rep_bwa = i
            break
    
    if rep_bwa == -1:
        return [BAD_SCORE] * 7
    
    # Parse CIGAR and compute scores
    seq_at = 0  # Position in read sequence
    ref_at = pos - REF_OFFSET  # Position in reference
    
    try:
        for length, op in parse_cigar(cigar):
            if op in ['S', 'I']:  # Soft clip or insertion
                seq_at += length
            elif op == 'D':  # Deletion
                ref_at += length
            elif op == 'M':  # Match/mismatch
                for i in range(length):
                    if seq_at >= len(seq) or ref_at < 0:
                        break
                    
                    qual_val = ord(qual[seq_at])
                    if qual_val >= MIN_QUAL:
                        # Quality score confidence calculation
                        if qual_val == ord(':') or qual_val == ord('F'):
                            conf = qual_val - ord('!')
                        else:
                            conf = (qual_val - ord('!')) // 2
                        
                        # Check mismatch against each repeat position
                        for rep in range(7):
                            ref_pos = ref_at + STARTS[rep] - STARTS[rep_bwa]
                            if 0 <= ref_pos < len(ref):
                                if seq[seq_at] != ref[ref_pos]:
                                    scores[rep] += conf
                            else:
                                scores[rep] = BAD_SCORE
                    
                    seq_at += 1
                    ref_at += 1
            else:
                # Unknown CIGAR operation
                return [BAD_SCORE] * 7
                
    except (ValueError, IndexError):
        return [BAD_SCORE] * 7
    
    return scores


def classify_variant(combined_scores):
    """
    Classify variant based on combined scores from read pair.
    
    Returns: '1B_KIV3', '1B_KIV2', '1B_tied', or '1A'
    """
    s = combined_scores
    
    # Check if position 0 (161032593) has minimum score -> 1B_KIV3
    if all(s[0] < s[i] for i in range(1, 7)):
        return '1B_KIV3'
    
    # Check if position 4 (161054784) has minimum score -> 1B_KIV2
    elif all(s[4] < s[i] for i in range(7) if i != 4):
        return '1B_KIV2'
    
    # Check if positions 0 and 4 are tied for minimum -> 1B_tied
    elif s[4] == s[0] and all(s[4] < s[i] for i in range(1, 7) if i != 4):
        return '1B_tied'
    
    # Otherwise check for 1A (another position has minimum)
    else:
        min_1A = BAD_SCORE
        for rep in range(1, 7):
            if rep != 4:
                min_1A = min(min_1A, s[rep])
        
        if min_1A < s[0] and min_1A < s[4]:
            return '1A'
    
    return None  # Ambiguous classification


def process_samtools_input(ref_seq):
    """
    Process samtools output from stdin.
    Expected format: qname pos cigar seq qual (5 columns, space-separated)
    """
    counts = {
        '1B_KIV3': 0,
        '1B_KIV2': 0,
        '1B_tied': 0,
        '1A': 0
    }
    
    prev_read = None
    prev_scores = None
    
    for line in sys.stdin:
        fields = line.strip().split()
        if len(fields) != 5:
            continue
        
        qname, pos_str, cigar, seq, qual = fields
        try:
            pos = int(pos_str)
        except ValueError:
            continue
        
        # Compute scores for this read
        scores = compute_scores(pos, cigar, seq, qual, ref_seq)
        
        # If this is the second read in a pair (same qname)
        if prev_read == qname:
            # Combine scores from both reads in pair
            combined = [min(BAD_SCORE, s1 + s2) for s1, s2 in zip(prev_scores, scores)]
            
            # Classify the variant
            variant_type = classify_variant(combined)
            if variant_type:
                counts[variant_type] += 1
            
            # Reset for next pair
            prev_read = None
            prev_scores = None
        else:
            # Store first read of pair
            prev_read = qname
            prev_scores = scores
    
    return counts


def main():
    parser = argparse.ArgumentParser(
        description='Realign LPA reads to classify variant types'
    )
    parser.add_argument(
        '--reference',
        required=True,
        help='Path to reference FASTA file (output from extract_reference)'
    )
    parser.add_argument(
        '--sample_id',
        required=True,
        help='Sample ID to prepend to output'
    )
    
    args = parser.parse_args()
    
    # Load reference sequence
    try:
        ref_seq = load_reference(args.reference)
    except FileNotFoundError:
        print(f"ERROR: Reference file not found: {args.reference}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"ERROR: Failed to load reference: {e}", file=sys.stderr)
        sys.exit(1)
    
    # Process samtools input from stdin
    counts = process_samtools_input(ref_seq)
    
    # Output results (tab-separated: sample_id, 1B_KIV3, 1B_KIV2, 1B_tied, 1A)
    print(f"{args.sample_id}\t{counts['1B_KIV3']}\t{counts['1B_KIV2']}\t{counts['1B_tied']}\t{counts['1A']}")


if __name__ == '__main__':
    main()