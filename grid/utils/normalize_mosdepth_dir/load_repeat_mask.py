# grid/utils/normalize_mosdepth_dir/load_repeat_mask.py
# In[1]: Imports
from collections import defaultdict

# In[2]: Function to Load Repeat Mask
def load_repeat_mask(repeat_bed: str) -> dict[str, set[int]]:
    """
    Load repeat regions into {chrom: set(kb_bins)}.
    
    Args:
        repeat_bed: Path to repeat mask BED file
    
    Returns:
        Dictionary mapping chromosomes to sets of excluded kb bins
    """
    excluded = defaultdict(set)
    with open(repeat_bed) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split()
            if len(parts) < 3:
                continue
            chrom, start_str, end_str = parts[0], parts[1], parts[2]
            try:
                start, end = int(start_str), int(end_str)
            except ValueError:
                continue
            for kb in range(start // 1000, end // 1000 + 1):
                excluded[chrom].add(kb)
    return excluded