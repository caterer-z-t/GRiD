# grid/utils/compute_dipcn_dir/write_dipcn_output.py
# In[1]: Imports
from typing import Dict
from pathlib import Path

# In[2]: Function to write diploid copy number output
def write_dipcn_output(results: Dict[str, float], output_file: str) -> None:
    """
    Write diploid copy numbers to output file.
    
    Output format:
    ID\tdipCN
    sample1\t2.345678
    sample2\t1.987654
    ...
    
    Args:
        results (Dict[str, float]): Dictionary mapping sample_id to diploid CN
        output_file (str): Output file path
        
    Returns:
        None
    """
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    with open(output_path, 'w') as f:
        f.write("ID\tdipCN\n")
        for sample_id, dip_cn in sorted(results.items()):
            f.write(f"{sample_id}\t{dip_cn:.6f}\n")