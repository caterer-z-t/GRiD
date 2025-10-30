# grid/utils/count_reads_dir/count_reads_in_region.py
# In[1]: Imports
from pathlib import Path
import pysam

# In[2]: Function to Count Reads in a Region
from ..ensure_crai import ensure_crai

def count_reads_in_region(
    cram_file: str,
    ref_fasta: str,
    chrom: str,
    start: int,
    end: int,
    proper_flags: set[int]
) -> int:
    """
    Count properly paired reads in a genomic region using pysam.
    Ensures CRAI index exists.
    """
    ensure_crai(cram_file, ref_fasta)

    count = 0
    with pysam.AlignmentFile(cram_file, "rc", reference_filename=ref_fasta) as bam:
        for read in bam.fetch(chrom, start, end):
            if read.flag in proper_flags:
                count += 1
    return count
