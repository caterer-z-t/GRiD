# grid/utils/helper_dir/find_all_cram_files.py
# In[1]: Imports
from pathlib import Path
from typing import List

# In[2]: Function to find all CRAM files in a directory
def find_cram_files(cram_dir: str) -> List[str]:
    """
    Return all CRAM files in a directory.
    
    Args:
        cram_dir (str): Path to the directory containing CRAM files.

    Returns:
        List[str]: List of paths to CRAM files.
    """
    cram_dir = Path(cram_dir).expanduser()
    return [str(f) for f in cram_dir.glob("*.cram")]
