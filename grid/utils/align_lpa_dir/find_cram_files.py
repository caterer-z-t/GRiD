# grid/utils/align_lpa_dir/find_cram_files.py
# In[1]: Imports
from pathlib import Path
from typing import List

# In[2]: Function to find CRAM/BAM files
def find_cram_files(cram_dir: str) -> List[Path]:
    """
    Find all CRAM and BAM files in the specified directory.
    
    Args:
        cram_dir (str): Directory to search for CRAM/BAM files
        
    Returns:
        List[Path]: List of Path objects for found CRAM/BAM files with valid indices
    """
    cram_dir_path = Path(cram_dir).expanduser().resolve()
    
    cram_files = []
    
    # Find CRAM files
    for cram_file in cram_dir_path.glob("*.cram"):
        # Check if index exists (.crai)
        if cram_file.with_suffix('.cram.crai').exists():
            cram_files.append(cram_file)
    
    # Find BAM files
    for bam_file in cram_dir_path.glob("*.bam"):
        # Check if index exists (.bai)
        if bam_file.with_suffix('.bam.bai').exists():
            cram_files.append(bam_file)
    
    # Sort for consistent ordering
    cram_files.sort()
    
    return cram_files