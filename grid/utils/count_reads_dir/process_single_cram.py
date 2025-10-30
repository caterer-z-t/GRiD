# grid/utils/count_reads_dir/process_single_cram.py
# In[1]: Imports
from pathlib import Path
from .count_reads_in_region import count_reads_in_region
from rich.console import Console

# In[2]: Setting up console for rich output
console = Console()

# In[3]: Function to process a single CRAM file
def process_single_cram(
    cram_path: str,
    ref_fasta: str,
    chrom: str,
    start: int,
    end: int,
    proper_flags: set[int]
) -> tuple[str, int | str]:
    """
    Process a single CRAM file and count reads in region.
    
    Args:
        cram_path: Path to CRAM file
        ref_fasta: Path to reference genome FASTA
        chrom: Chromosome name
        start: Start position
        end: End position
        proper_flags: Set of SAM flags to count
    
    Returns:
        Tuple of (basename, count) where count is int or "Error"
    """
    basename = Path(cram_path).name
    
    try:
        count = count_reads_in_region(
            cram_path, 
            ref_fasta, 
            chrom, 
            start, 
            end, 
            proper_flags
        )
        return (basename, count)
    except Exception as e:
        console.print(f"[red]✗ Failed to count reads for {basename}: {e}[/red]")
        return (basename, "Error")