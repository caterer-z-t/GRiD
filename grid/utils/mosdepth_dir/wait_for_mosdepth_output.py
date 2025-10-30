# grid/utils/mosdepth_dir/wait_for_mosdepth_output.py
# In[1]: Imports
from pathlib import Path
import time
from rich.console import Console

# In[2]: Console Setup
console = Console()

# In[3]: Function to Wait for Mosdepth Output
def wait_for_mosdepth_output(
    work_dir: Path,
    sample_name: str,
    max_attempts: int = 3,
    sleep_seconds: int = 2
) -> Path:
    """
    Wait for mosdepth output file to appear.
    
    Args:
        work_dir: Working directory where mosdepth writes output
        sample_name: Sample name to search for
        max_attempts: Maximum number of attempts to find file
        sleep_seconds: Seconds to wait between attempts
    
    Returns:
        Path to regions.bed.gz file
    
    Raises:
        FileNotFoundError: If file is not found after max attempts
    """
    for attempt in range(max_attempts):
        matches = list(work_dir.glob(f"{sample_name}*regions.bed.gz"))
        if matches:
            return matches[0]
        
        console.print(f"[yellow]Waiting for mosdepth output ({attempt + 1}/{max_attempts})...[/yellow]")
        time.sleep(sleep_seconds)
    
    raise FileNotFoundError(
        f"Expected mosdepth output matching pattern '{sample_name}*regions.bed.gz' "
        f"not found in {work_dir} for sample {sample_name}."
    )