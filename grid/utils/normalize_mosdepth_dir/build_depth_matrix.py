# grid/utils/normalize_mosdepth_dir/build_depth_matrix.py
# In[1]: Imports
from collections import defaultdict
import numpy as np

# In[2]: Function to Build Depth Matrix
def build_depth_matrix(
    regions_to_extract: dict[str, list[tuple[int, int, float]]]
) -> dict[str, dict[tuple[int, int], float]]:
    """
    Build depth matrix from extracted regions.
    
    Args:
        regions_to_extract: Dictionary of sample IDs to region lists
    
    Returns:
        Dictionary mapping {sample: {(start,end): depth}}
    """
    depth_matrix = defaultdict(dict)
    for ind, regions in regions_to_extract.items():
        for start, end, depth in regions:
            depth_matrix[ind][(start, end)] = depth
    return depth_matrix

# In[3]: Function to Build Matrix Arrays
def build_matrix_from_regions(regions_to_extract, individuals_order=None):
    """
    Build numpy matrix from regions_to_extract.

    Args:
        regions_to_extract: dict of {individual_id: [(start,end,depth),...]}
        individuals_order: optional list of individual IDs to order rows
    
    Returns:
        individuals_order: list of individual IDs in order
        regions_list: list of (start,end) tuples in order
        mat: numpy array of shape (n_individuals, n_regions) with depths
    """
    if individuals_order is None:
        individuals_order = sorted(regions_to_extract.keys())
    # Collect all regions as tuples (start,end)
    region_set = set()
    for ind in individuals_order:
        for s, e, d in regions_to_extract.get(ind, []):
            region_set.add((s, e))
    regions_list = sorted(region_set)

    n_inds = len(individuals_order)
    n_regs = len(regions_list)
    mat = np.full((n_inds, n_regs), np.nan, dtype=float)

    region_index = {r: i for i, r in enumerate(regions_list)}
    ind_index = {ind: i for i, ind in enumerate(individuals_order)}

    # Fill matrix with raw depths
    for ind, regions in regions_to_extract.items():
        i = ind_index[ind]
        for s, e, d in regions:
            j = region_index.get((s, e))
            if j is not None:
                mat[i, j] = d

    return individuals_order, regions_list, mat