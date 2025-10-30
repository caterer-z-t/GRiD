# grid/utils/mosdepth_dir/build_mosdepth_command.py
# In[1]: Imports
from pathlib import Path
from typing import List

# In[2]: Function to build mosdepth command
def build_mosdepth_command(
    cram_path: str,
    ref_fasta: str,
    output_prefix: Path,
    by: int,
    fast_mode: bool
) -> List[str]:
    """
    Build mosdepth command with appropriate flags.
    
    Args:
        cram_path: Path to CRAM file
        ref_fasta: Path to reference genome FASTA
        output_prefix: Output prefix for mosdepth
        by: Bin size for coverage calculation
        fast_mode: Whether to use fast mode
    
    Returns:
        List of command arguments
    """
    cmd = ["mosdepth", "-n", "--by", str(by), "-f", str(ref_fasta), str(output_prefix), str(cram_path)]
    
    if fast_mode:
        cmd.insert(1, "--fast-mode")
    
    return cmd