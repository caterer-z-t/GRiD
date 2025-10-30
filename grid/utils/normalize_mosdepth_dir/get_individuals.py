# grid/utils/normalize_mosdepth_dir/get_individuals.py
# In[1]: Imports
from pathlib import Path

# In[2]: Function to Get Individuals
def get_individuals(mosdepth_dir: str) -> dict[str, Path]:
    """
    Map sample IDs to their mosdepth.global.dist files.
    
    Args:
        mosdepth_dir: Directory containing mosdepth output files
    
    Returns:
        Dictionary mapping sample IDs to their global.dist file paths
    """
    mosdepth_dir = Path(mosdepth_dir)
    return {
        f.stem.split(".")[0]: f
        for f in mosdepth_dir.glob("*.mosdepth.global.dist.txt")
    }