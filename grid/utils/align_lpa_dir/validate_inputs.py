# grid/utils/align_lpa_dir/validate_inputs.py
# In[1]: Imports
from pathlib import Path

# In[2]: Function to validate input files and directories
def validate_inputs(cram_dir, reference_fa, lpa_ref_fasta, positions_file):
    """
    Validate that all required input files and directories exist.
    
    Args:
        cram_dir (str): Directory containing CRAM/BAM files
        reference_fa (str): Reference genome FASTA
        lpa_ref_fasta (str): LPA reference FASTA
        positions_file (str): Positions file
        
    Raises:
        FileNotFoundError: If any required file/directory is missing
    """
    # Check CRAM directory
    cram_path = Path(cram_dir).expanduser().resolve()
    if not cram_path.exists():
        raise FileNotFoundError(f"CRAM directory not found: {cram_dir}")
    if not cram_path.is_dir():
        raise NotADirectoryError(f"CRAM path is not a directory: {cram_dir}")
    
    # Check reference FASTA
    ref_path = Path(reference_fa).expanduser().resolve()
    if not ref_path.exists():
        raise FileNotFoundError(f"Reference FASTA not found: {reference_fa}")
    
    # Check LPA reference FASTA
    lpa_ref_path = Path(lpa_ref_fasta).expanduser().resolve()
    if not lpa_ref_path.exists():
        raise FileNotFoundError(f"LPA reference FASTA not found: {lpa_ref_fasta}")
    
    # Check positions file
    pos_path = Path(positions_file).expanduser().resolve()
    if not pos_path.exists():
        raise FileNotFoundError(f"Positions file not found: {positions_file}")