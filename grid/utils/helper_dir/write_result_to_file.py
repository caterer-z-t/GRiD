# grid/utils/helper_dir/write_result_to_file.py
# In[1]: Imports
from pathlib import Path
from threading import Lock
from rich.console import Console

# In[2]: Initialize console for rich output
console = Console()

# In[3]: Define function to write results to file
def write_result_to_file(
    output_file: Path,
    basename: str,
    count: int | str,
    write_lock: Lock,
    progress_console=None
) -> None:
    """
    Write a single result to the output file in a thread-safe manner.
    
    Args:
        output_file: Path to output file
        basename: Sample basename
        count: Read count or "Error"
        write_lock: Threading lock for safe file writing
        progress_console: Optional Progress console for proper rendering
    """
    with write_lock:
        with open(output_file, "a") as f:
            f.write(f"{basename}\t{count}\n")
    
    # Display result after writing - use progress console if available
    if count != "Error":
        if progress_console:
            progress_console.print(f"[green]✓ {basename} done: reads={count}[/green]")
        else:
            console.print(f"[green]✓ {basename} done: reads={count}[/green]")