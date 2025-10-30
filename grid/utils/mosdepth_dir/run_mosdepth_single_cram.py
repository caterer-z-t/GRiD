# grid/utils/mosdepth_dir/run_mosdepth_single_cram.py
# In[1]: Imports
from pathlib import Path
import subprocess
from typing import Tuple

from ..ensure_crai import ensure_crai
from .build_mosdepth_command import build_mosdepth_command
from .wait_for_mosdepth_output import wait_for_mosdepth_output
from .compute_region_coverage import compute_region_coverage

# In[2]: Function to Run Mosdepth on Single CRAM
def run_mosdepth_single_cram(
    cram_path: str,
    ref_fasta: str,
    work_dir: Path,
    chrom: str,
    start: int,
    end: int,
    region_name: str,
    by: int,
    fast_mode: bool
) -> Tuple[str, int | str]:
    """
    Run mosdepth on a single CRAM file and compute coverage for region.
    
    Args:
        cram_path: Path to CRAM file
        ref_fasta: Path to reference genome FASTA
        work_dir: Working directory for intermediate files
        chrom: Chromosome name
        start: Start position
        end: End position
        by: Bin size for mosdepth
        fast_mode: Whether to use fast mode
    
    Returns:
        Tuple of (sample_name, coverage) where coverage is int or "Error"
    """
    cram = Path(cram_path)
    sample_name = cram.stem
    
    try:
        # Ensure CRAI index exists
        ensure_crai(str(cram), ref_fasta)
        
        # Build and run mosdepth command
        out_prefix = work_dir / f"{sample_name}_{region_name}"
        cmd = build_mosdepth_command(str(cram), ref_fasta, out_prefix, by, fast_mode)
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        
        # Wait for output file
        regions_file = wait_for_mosdepth_output(work_dir, sample_name)
        
        # Compute coverage
        coverage = compute_region_coverage(regions_file, chrom, start, end)
        
        return (sample_name, coverage)
        
    except Exception as e:
        return (sample_name, "Error")