# grid/utils/align_lpa_dir/load_reference.py
# In[1]: Imports
from pathlib import Path

# In[2]: Function to load reference sequence from FASTA
def load_reference(ref_fasta_path: str) -> str:
    """
    Load reference sequence from FASTA file.
    
    Args:
        ref_fasta_path (str): Path to reference FASTA file
        
    Returns:
        str: Reference sequence (without header lines)
        
    Raises:
        FileNotFoundError: If FASTA file doesn't exist
        ValueError: If FASTA file is empty or malformed
    """
    fasta_path = Path(ref_fasta_path).expanduser().resolve()
    
    with open(fasta_path, 'r') as f:
        lines = f.readlines()
    
    if len(lines) == 0:
        raise ValueError(f"Empty FASTA file: {ref_fasta_path}")
    
    # Skip header line(s) and concatenate sequence
    ref_seq = ''.join(line.strip() for line in lines if not line.startswith('>'))
    
    if len(ref_seq) == 0:
        raise ValueError(f"No sequence found in FASTA file: {ref_fasta_path}")
    
    return ref_seq