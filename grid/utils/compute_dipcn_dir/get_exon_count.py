# grid/utils/compute_dipcn_dir/get_exon_count.py
# In[1]: Imports
from typing import Dict

# In[2]: Function to get exon count
def get_exon_count(counts_dict: Dict[str, int], exon_type: str) -> int:
    """
    Get the appropriate count for a given exon type.
    
    Exon types:
    - 1B_KIV3: just the 1B_KIV3 count
    - 1B_notKIV3: 1B_KIV2 + 1B_tied
    - 1B: 1B_KIV3 + 1B_KIV2 + 1B_tied
    - 1A: just the 1A count
    
    Args:
        counts_dict (Dict[str, int]): Dictionary with keys '1B_KIV3', '1B_KIV2', '1B_tied', '1A'
        exon_type (str): Exon type to compute
        
    Returns:
        int: Count for the specified exon type
        
    Raises:
        ValueError: If exon_type is not recognized
    """
    if exon_type == '1B_KIV3':
        return counts_dict.get('1B_KIV3', 0)
    
    elif exon_type == '1B_notKIV3':
        return counts_dict.get('1B_KIV2', 0) + counts_dict.get('1B_tied', 0)
    
    elif exon_type == '1B':
        return (
            counts_dict.get('1B_KIV3', 0) + 
            counts_dict.get('1B_KIV2', 0) + 
            counts_dict.get('1B_tied', 0)
        )
    
    elif exon_type == '1A':
        return counts_dict.get('1A', 0)
    
    else:
        raise ValueError(f"Unknown exon type: {exon_type}")