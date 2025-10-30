# grid/utils/normalize_mosdepth_dir/select_high_variance_regions.py
# In[1]: Function to Select High Variance Regions
def select_high_variance_regions(variance_ratios: dict[int, float], top_frac: float = 0.1) -> list[int]:
    """
    Select top fraction of regions by variance ratio.

    Args:
        variance_ratios: {region_index: variance_ratio}
        top_frac: fraction of top regions to retain

    Returns:
        List of selected region indices
    """
    if not variance_ratios:
        return []
    sorted_indices = sorted(variance_ratios.items(), key=lambda x: x[1], reverse=True)
    n_keep = max(1, int(len(sorted_indices) * top_frac))
    return [idx for idx, _ in sorted_indices[:n_keep]]
