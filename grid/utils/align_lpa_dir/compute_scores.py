# grid/utils/align_lpa_dir/compute_scores.py
# In[1]: Imports
import re
from typing import List, Tuple

# In[2]: Constants and Helper Functions
# Constants
MIN_QUAL = ord(':')  # ASCII value for ':'
BAD_SCORE = 9999

# In[3]: Function to parse CIGAR string
def parse_cigar(cigar: str) -> List[Tuple[int, str]]:
    """
    Parse CIGAR string into list of (length, operation) tuples.
    
    Args:
        cigar (str): CIGAR string (e.g., '100M', '50M1I49M')
        
    Returns:
        List[Tuple[int, str]]: List of (length, operation) tuples
    """
    pattern = re.compile(r'(\d+)([MIDNSHP=X])')
    return [(int(length), op) for length, op in pattern.findall(cigar)]

# In[4]: Function to find closest repeat
def find_closest_repeat(pos: int, starts: List[int], window: int = 1000) -> int:
    """
    Find the closest repeat to a given position within a window.
    
    Args:
        pos (int): Read position
        starts (List[int]): List of repeat start positions
        window (int): Maximum distance to consider (default: 1000bp)
        
    Returns:
        int: Index of closest repeat, or -1 if none within window
    """
    min_dist = float('inf')
    closest_rep = -1
    
    for i, start in enumerate(starts):
        dist = abs(pos - start)
        if dist < min_dist:
            min_dist = dist
            closest_rep = i
    
    # Return -1 if outside window
    if min_dist > window:
        return -1
    
    return closest_rep

def compute_alignment_scores(
    pos: int,
    cigar: str,
    seq: str,
    qual: str,
    ref_seq: str,
    starts: List[int],
    ref_offset: int
) -> List[int]:
    """
    Compute mismatch scores for each of the 7 repeat positions.
    
    Lower score = better match to that repeat position.
    
    UPDATED: Now uses ±1000bp window and closest repeat matching.
    
    Args:
        pos (int): Read alignment position (1-based from BAM)
        cigar (str): CIGAR string
        seq (str): Read sequence
        qual (str): Quality string
        ref_seq (str): Reference sequence
        starts (List[int]): List of 7 repeat start positions (1-based)
        ref_offset (int): Reference offset position (1-based)
        
    Returns:
        List[int]: Array of 7 scores, one for each potential repeat alignment
    """
    scores = [0] * 7
    
    # Find which repeat this read likely aligns to based on position
    # UPDATED: Use closest repeat within 1000bp instead of strict ±500bp
    rep_bwa = find_closest_repeat(pos, starts, window=1000)
    
    if rep_bwa == -1:
        return [BAD_SCORE] * 7
    
    # Parse CIGAR and compute scores
    seq_at = 0  # Position in read sequence
    ref_at = pos - ref_offset  # Position in reference
    
    try:
        for length, op in parse_cigar(cigar):
            if op in ['S', 'H']:  # Soft clip or hard clip
                if op == 'S':
                    seq_at += length
                # Hard clips don't consume query
            elif op == 'I':  # Insertion
                seq_at += length
            elif op == 'D':  # Deletion
                ref_at += length
            elif op == 'N':  # Skipped region
                ref_at += length
            elif op in ['M', '=', 'X']:  # Match/mismatch
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
                            ref_pos = ref_at + starts[rep] - starts[rep_bwa]
                            if 0 <= ref_pos < len(ref_seq):
                                if seq[seq_at] != ref_seq[ref_pos]:
                                    scores[rep] += conf
                            else:
                                scores[rep] = BAD_SCORE
                    
                    seq_at += 1
                    ref_at += 1
            else:
                # Unknown CIGAR operation - should not happen
                pass
                
    except (ValueError, IndexError):
        return [BAD_SCORE] * 7
    
    return scores