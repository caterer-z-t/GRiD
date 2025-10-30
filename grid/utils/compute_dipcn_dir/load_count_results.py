# grid/utils/compute_dipcn_dir/load_count_results.py
# In[1]: Imports
from pathlib import Path
from typing import Dict
from .normalize_sample_id import normalize_sample_id

# In[2]: Function to load count results
def load_count_results(count_file: Path) -> Dict[str, Dict[str, int]]:
    """
    Load read counts from realignment output.
    
    Expected format: sample_id\t1B_KIV3\t1B_KIV2\t1B_tied\t1A
    
    Args:
        count_file (Path): Path to count file
        
    Returns:
        Dict[str, Dict[str, int]]: Dictionary mapping sample_id to count dict
            Example: {'NWD123': {'1B_KIV3': 45, '1B_KIV2': 23, ...}}
    """
    counts = {}
    
    with open(count_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
                
            fields = line.split('\t')
            if len(fields) != 5:
                continue
            
            # Normalize sample ID
            original_id = fields[0]
            sample_id = normalize_sample_id(original_id)
            
            try:
                counts[sample_id] = {
                    '1B_KIV3': int(fields[1]),
                    '1B_KIV2': int(fields[2]),
                    '1B_tied': int(fields[3]),
                    '1A': int(fields[4])
                }
            except ValueError:
                # Skip malformed lines
                continue
    
    return counts