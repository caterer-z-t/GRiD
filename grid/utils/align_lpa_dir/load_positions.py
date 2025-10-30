# grid/utils/align_lpa_dir/load_positions.py
# In[1]: Imports
from pathlib import Path
from typing import Tuple, List

# In[2]: Function to load positions from file
def load_positions(positions_file: str, genome_build: str) -> Tuple[List[int], int]:
    """
    Load hardcoded positions from file.
    
    File format (tab-separated):
    Hg38file    Hg19file
    160611000   161032032
    160611561   161032593
    ...
    
    Args:
        positions_file (str): Path to positions file
        genome_build (str): Genome build ('hg38' or 'hg19')
        
    Returns:
        Tuple[List[int], int]: (list of 7 repeat starts, reference offset)
        
    Raises:
        ValueError: If positions file is malformed or doesn't contain 7 repeats
    """
    positions_path = Path(positions_file).expanduser().resolve()
    
    starts = []
    ref_offset = 0
    
    # Determine which column to use (0=hg38, 1=hg19)
    col_idx = 0 if genome_build == 'hg38' else 1
    
    with open(positions_path, 'r') as f:
        # Skip header
        header = next(f)
        
        for i, line in enumerate(f):
            line = line.strip()
            if not line:
                continue
                
            fields = line.split('\t')
            if len(fields) < 2:
                raise ValueError(f"Invalid line in positions file: {line}")
            
            try:
                pos = int(fields[col_idx])
                
                # First position is the reference offset
                if i == 0:
                    ref_offset = pos
                # Remaining positions are the repeat starts
                else:
                    starts.append(pos)
            except ValueError as e:
                raise ValueError(f"Could not parse position from line: {line}") from e
    
    if len(starts) != 7:
        raise ValueError(f"Expected 7 repeat positions, got {len(starts)}")
    
    return starts, ref_offset