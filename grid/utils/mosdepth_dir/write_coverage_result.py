# grid/utils/mosdepth_dirwrite_coverage_result.py
# In[1]: Imports
from pathlib import Path
from threading import Lock
import csv
from ..helper_dir.display_results import print_individual_success

# In[2]: Function to Write Coverage Result
def write_coverage_result(
    output_file: Path,
    sample_name: str,
    coverage: int,
    write_lock: Lock,
    progress_console=None
) -> None:
    """
    Write coverage result to output file in a thread-safe manner.
    
    Args:
        output_file: Path to output TSV file
        sample_name: Sample name
        coverage: Coverage value
        write_lock: Threading lock for safe file writing
        progress_console: Optional Progress console for proper rendering
    """
    with write_lock:
        with open(output_file, "a", newline="") as f:
            writer = csv.writer(f, delimiter="\t")
            writer.writerow([sample_name, coverage])
    
    # Display result
    print_individual_success(sample_name, f"done: coverage={coverage}", progress_console)