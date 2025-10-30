# grid/utils/mosdepth_dir/compute_region_coverage.py
# In[1]: Imports
import gzip
from pathlib import Path

# In[2]: Function to Compute Region Coverage
def compute_region_coverage(
    regions_file: Path,
    chrom: str,
    start: int,
    end: int
) -> int:
    """
    Compute average coverage for a genomic region from mosdepth output.
    
    Args:
        regions_file: Path to mosdepth regions.bed.gz file
        chrom: Chromosome name
        start: Start position
        end: End position
    
    Returns:
        Average coverage as integer (rounded and scaled by 100)
    """
    region_cov = 0.0
    covered_bp = 0

    with gzip.open(regions_file, "rt") as f:
        for line in f:
            r_chr, r_start, r_end, mean_cov = line.strip().split("\t")[:4]
            r_start, r_end = int(r_start), int(r_end)
            mean_cov = float(mean_cov)

            if r_chr != chrom:
                continue

            overlap_start = max(start, r_start)
            overlap_end = min(end, r_end)
            overlap = overlap_end - overlap_start

            if overlap > 0:
                region_cov += mean_cov * overlap
                covered_bp += overlap

    return int(round(100 * (region_cov / covered_bp))) if covered_bp > 0 else 0