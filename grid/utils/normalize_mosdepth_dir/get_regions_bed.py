# grid/utils/normalize_mosdepth_dir/get_regions_bed.py
# In[1]: Imports
import gzip
from pathlib import Path

# In[2]: Function to Get Regions BED Path
def regions_bed_gz(global_dist_file: Path) -> Path:
    """
    Get the path to the regions BED file corresponding to a global distribution file.

    Args:
        global_dist_file: Path to the mosdepth.global.dist.txt file

    Returns:
        Path to the corresponding .regions.bed.gz file
    """
    return global_dist_file.with_name(
        global_dist_file.name.replace(".mosdepth.global.dist.txt", ".regions.bed.gz")
    )