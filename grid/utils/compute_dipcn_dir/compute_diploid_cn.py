# grid/utils/compute_dipcn_dir/compute_diploid_cn.py
# In[1]: Imports
from typing import Dict, List, Tuple
from .get_exon_count import get_exon_count

# In[2]: Function to compute diploid copy number
def compute_diploid_cn_for_exon(
    counts: Dict[str, Dict[str, int]],
    neighbors: Dict[str, Tuple[float, List[Tuple[str, float, float]]]],
    exon_type: str,
    n_neighbors: int = 200
) -> Dict[str, float]:
    """
    Compute diploid copy number for all samples for a given exon type.
    
    Formula: dipCN = (sample_count / sample_scale) / (mean_neighbor_normalized_count)
    where mean_neighbor_normalized_count = mean(neighbor_count / neighbor_scale)
    
    Args:
        counts (Dict): Dictionary mapping sample_id to count dict
        neighbors (Dict): Dictionary mapping sample_id to (scale, neighbor_list)
        exon_type (str): Exon type to compute ('1B_KIV3', '1B_notKIV3', '1B', '1A')
        n_neighbors (int): Number of top neighbors to use
        
    Returns:
        Dict[str, float]: Dictionary mapping sample_id to diploid copy number
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
        
        for neighbor_id, neighbor_scale, distance in neighbor_list[:n_neighbors]:
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