# grid/utils/helper_dir/create_region.py
# In[1]: Imports
import sys
from rich.console import Console

console = Console()

# In[2]: Define function to create region string
def create_region_string(region: str, chrom: str, start: int, end: int) -> str:
    """
    Create a genomic region string in 'chr:start-end' format.

    Args:
        region: Full region string (e.g., 'chr6:160000000-160100000')
        chrom: Chromosome name (e.g., 'chr6')
        start: Start position (e.g., 160000000)
        end: End position (e.g., 160100000)

    Returns:
        Genomic region string.
    """
    if not region:
        if chrom and start is not None and end is not None:
            region = f"{chrom}:{start}-{end}"
        else:
            console.print("[red]✗ You must provide either a full region or chromosome, start, and end separately.[/red]")
            sys.exit(1)
    return region