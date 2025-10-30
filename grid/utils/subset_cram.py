# grid/utils/subset_cram.py
# In[1]: Imports
from pathlib import Path
import pysam

# In[2]: Define function to subset CRAM file
def subset_cram(cram_path: str, region: str, output_path: str, reference: str = None) -> str:
    """
    Subset a CRAM file for a specific genomic region.

    Args:
        cram_path: Input CRAM file.
        region: Region in 'chr:start-end' format (e.g., 'chr6:160000000-160100000').
        output_path: Path to output CRAM file.
        reference: Reference genome FASTA (required for CRAM output)

    Returns:
        Path to the subset CRAM file.
    """
    cram_file = Path(cram_path)
    if not cram_file.exists():
        raise FileNotFoundError(f"CRAM file not found: {cram_path}")

    out_file = Path(output_path)
    with pysam.AlignmentFile(str(cram_file), "rc", reference_filename=reference) as in_bam, \
         pysam.AlignmentFile(str(out_file), "wc", template=in_bam, reference_filename=reference) as out_bam:
        for read in in_bam.fetch(region=region):
            out_bam.write(read)

    return str(out_file)
