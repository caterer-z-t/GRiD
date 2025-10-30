# grid/utils/compute_dipcn_dir/validate_sample_overlap.py
# In[1]: Imports
from typing import Dict, Set
from rich.console import Console

# In[2]: Function to validate sample overlap
def validate_sample_overlap(
    counts: Dict,
    neighbors: Dict,
    console: Console = None
) -> int:
    """
    Validate that there are overlapping samples between counts and neighbors.
    
    Args:
        counts (Dict): Dictionary of count results
        neighbors (Dict): Dictionary of neighbor results
        console (Console): Rich console for output (optional)
        
    Returns:
        int: Number of overlapping samples
    """
    count_samples = set(counts.keys())
    neighbor_samples = set(neighbors.keys())
    
    overlap = count_samples & neighbor_samples
    
    if console:
        console.print(f"  • Samples in count file: {len(count_samples)}")
        console.print(f"  • Samples in neighbor file: {len(neighbor_samples)}")
        console.print(f"  • Overlapping samples: {len(overlap)}")
    
    return len(overlap)