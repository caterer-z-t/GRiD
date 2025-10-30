# grid/utils/helper_dir/setup_output_file.py
# In[1]: Imports
from pathlib import Path

# In[2]: Function to set up output file
def setup_output_file(output_file: str, chrom: str, start: int, end: int) -> Path:
    """
    Create output file with header.
    
    Args:
        output_file: Path to output TSV file
        chrom: Chromosome name
        start: Start position
        end: End position
    
    Returns:
        Path object for the output file
    """
    output_path = Path(output_file).expanduser()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    with open(output_path, "w") as f:
        f.write(f"Sample\t{chrom}:{start}-{end}\n")
    
    return output_path