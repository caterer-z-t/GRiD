# grid/utils/find_neighbors_dir/read_normalized_data.py
# In[1]: Imports
import gzip
import numpy as np
from rich.console import Console

# In[2]: Console setup
console = Console()

# In[3]: Function to read normalized data
def read_normalized_data(input_file):
    """
    Read gzipped normalized depth data.
    
    Args:
        input_file (str): Path to gzipped normalized depth file.

    Returns:
        individuals (list): List of individual IDs.
        regions (list): List of region identifiers.
        mat (np.ndarray): Data matrix of shape [regions x individuals].
        scales (dict): Dictionary of scaling factors for each individual.
    """
    console.print(f"[yellow]Reading normalized data from {input_file}[/yellow]")

    with gzip.open(input_file, 'rt') as f:
        header = f.readline().strip().split('\t')
        individuals = header[1:]
        regions, data = [], []

        for line in f:
            parts = line.strip().split('\t')
            regions.append(parts[0])
            vals = [
                float(v) if v not in ('NA', 'nan') else np.nan
                for v in parts[1:]
            ]
            data.append(vals)

    mat = np.array(data)
    scales = {ind: 1.0 for ind in individuals}

    console.print(f"[green]✓ Loaded matrix of shape {mat.shape}[/green]")
    return individuals, regions, mat, scales
