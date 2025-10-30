# grid/utils/compute_dipcn_dir/load_neighbor_results.py
# In[1]: Imports
import gzip
from pathlib import Path
from typing import Dict, List, Tuple
from .normalize_sample_id import normalize_sample_id

# In[2]: Function to load neighbor results
def load_neighbor_results(neighbor_file: Path) -> Dict[str, Tuple[float, List[Tuple[str, float, float]]]]:
    """
    Load neighbor information from find_neighbors output.
    
    Expected format: sample_id\tscale\tneighbor1\tneighbor1_scale\tneighbor1_distance\t...
    
    Args:
        neighbor_file (Path): Path to neighbor file (can be gzipped)
        
    Returns:
        Dict mapping sample_id to (scale, neighbor_list)
        where neighbor_list is [(neighbor_id, neighbor_scale, distance), ...]
        Example: {'NWD123': (1.0, [('NWD456', 0.98, 0.05), ('NWD789', 1.02, 0.07), ...])}
    """
    neighbors = {}
    
    # Handle both gzipped and plain text files
    if str(neighbor_file).endswith('.gz'):
        open_func = gzip.open
        mode = 'rt'
    else:
        open_func = open
        mode = 'r'
    
    with open_func(neighbor_file, mode) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
                
            fields = line.split('\t')
            if len(fields) < 2:
                continue
            
            # Normalize sample ID
            original_id = fields[0]
            sample_id = normalize_sample_id(original_id)
            
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