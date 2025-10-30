# grid/utils/align_lpa_dir/classify_reads.py
# In[1]: Imports
from typing import List, Dict
from collections import defaultdict
from .compute_scores import compute_alignment_scores
from .classify_variant import classify_variant_from_scores

# In[2]: Function to classify read pairs
def classify_read_pairs(
    reads: List[Dict],
    ref_seq: str,
    starts: List[int],
    ref_offset: int
) -> Dict[str, int]:
    """
    Classify read pairs into variant categories.
    
    UPDATED: Now properly groups reads by qname before processing pairs.
    This fixes the bug where pairs had to be consecutive in the list.
    
    Args:
        reads (List[Dict]): List of read dictionaries
        ref_seq (str): Reference sequence
        starts (List[int]): List of repeat start positions
        ref_offset (int): Reference offset
        
    Returns:
        Dict[str, int]: Counts for each variant type:
            - '1B_KIV3': Count of 1B_KIV3 variants
            - '1B_KIV2': Count of 1B_KIV2 variants
            - '1B_tied': Count of 1B_tied variants
            - '1A': Count of 1A variants
    """
    counts = {
        '1B_KIV3': 0,
        '1B_KIV2': 0,
        '1B_tied': 0,
        '1A': 0
    }
    
    # Group reads by qname (read pair name)
    read_pairs = defaultdict(list)
    for read in reads:
        read_pairs[read['qname']].append(read)
    
    BAD_SCORE = 9999
    
    # Process each read pair
    for qname, pair_reads in read_pairs.items():
        # Only process proper pairs (exactly 2 reads)
        if len(pair_reads) != 2:
            continue
        
        # Compute scores for both reads in the pair
        scores_1 = compute_alignment_scores(
            pos=pair_reads[0]['pos'],
            cigar=pair_reads[0]['cigar'],
            seq=pair_reads[0]['seq'],
            qual=pair_reads[0]['qual'],
            ref_seq=ref_seq,
            starts=starts,
            ref_offset=ref_offset
        )
        
        scores_2 = compute_alignment_scores(
            pos=pair_reads[1]['pos'],
            cigar=pair_reads[1]['cigar'],
            seq=pair_reads[1]['seq'],
            qual=pair_reads[1]['qual'],
            ref_seq=ref_seq,
            starts=starts,
            ref_offset=ref_offset
        )
        
        # Skip if either read failed to score
        if scores_1[0] == BAD_SCORE or scores_2[0] == BAD_SCORE:
            continue
        
        # Combine scores from both reads in pair
        combined = [
            min(BAD_SCORE, s1 + s2) 
            for s1, s2 in zip(scores_1, scores_2)
        ]
        
        # Classify the variant
        variant_type = classify_variant_from_scores(combined)
        if variant_type:
            counts[variant_type] += 1
    
    return counts