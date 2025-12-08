# grid/utils/create_crai.py
# In[1]: Imports
from pathlib import Path
import pysam

# In[2]: Function to Ensure CRAI Index
def ensure_crai(cram_path: str, reference: str = None) -> str:
    """
    Ensure a CRAI index exists for the given CRAM file.

    Args:
        cram_path: Path to CRAM file.
        reference: Optional reference genome FASTA (required if CRAM is unindexed)

    Returns:
        Path to CRAI file.
    """
    cram_file = Path(cram_path)

    if not cram_file.exists():
        raise FileNotFoundError(f"CRAM not found: {cram_path}")

    # pysam creates file.crai, not file.cram.crai
    crai_file = cram_file.with_suffix(".crai")

    if not crai_file.exists():
        if reference is None:
            raise ValueError("Reference genome must be provided")
        pysam.index(str(cram_file), reference=reference)

    return cram_file
