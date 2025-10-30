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
        raise FileNotFoundError(f"CRAM file not found: {cram_path}")

    crai_file = cram_file.with_suffix(cram_file.suffix + ".crai")
    if crai_file.exists():
        return str(crai_file)

    if reference is None:
        raise ValueError("Reference genome must be provided to create CRAI.")

    pysam.index(str(cram_file), reference=reference)
    return str(crai_file)
