# grid/utils/align_lpa_dir/classify_variant.py
# In[1]: Imports
from typing import List, Optional

BAD_SCORE = 9999

# In[2]: Function to classify variant from scores
def classify_variant_from_scores(combined_scores: List[int]) -> Optional[str]:
    """
    Classify variant based on combined scores from read pair.
    
    Args:
        combined_scores (List[int]): Array of 7 combined scores
        
    Returns:
        Optional[str]: Variant type ('1B_KIV3', '1B_KIV2', '1B_tied', '1A') or None
    """
    s = combined_scores
    
    # Check if position 0 (first repeat) has minimum score -> 1B_KIV3
    if all(s[0] < s[i] for i in range(1, 7)):
        return '1B_KIV3'
    
    # Check if position 4 (fifth repeat) has minimum score -> 1B_KIV2
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