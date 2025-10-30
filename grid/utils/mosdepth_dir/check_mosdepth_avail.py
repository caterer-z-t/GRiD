# grid/utils/mosdepth_dir/check_mosdepth_avail.py
# In[1]: Imports
import shutil

# In[2]: Check Mosdepth Availability
def check_mosdepth_available():
    """
    Check if mosdepth executable is available in PATH.
    
    Args:
        None

    Returns:
        None
        
    Raises:
        RuntimeError: If mosdepth is not found.
    """
    if shutil.which("mosdepth") is None:
        raise RuntimeError("mosdepth not found in PATH. Please install mosdepth before running.")