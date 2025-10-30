# grid/utils/find_neighbors_dir/save_neighbors.py
# In[1]: Imports
import gzip
from rich.console import Console

# In[2]: Set up console
console = Console()

# In[3]: Function to save neighbors
def save_neighbors(neighbors_dict, scales, output_prefix, zmax, R_use=None):
    """
    Save the nearest neighbors dictionary to a gzipped text file.
    
    Args:
        neighbors_dict: dict of {individual_id: [(neighbor_id, distance), ...]}
        scales: dict of {individual_id: scale_value}
        output_prefix: str, prefix for output file
        zmax: float, max z-score clipping
        R_use: int, number of regions used for normalization (for normalized distance)
    """
    if R_use is None:
        # If not provided, assume 1 to avoid division by zero
        R_use = 1

    output_file = f"{output_prefix}.zMax{zmax:.1f}.txt.gz"
    with gzip.open(output_file, "wt") as out:
        for ind, neighbors in neighbors_dict.items():
            out.write(f"{ind}\t{scales.get(ind, 1.0):.2f}")
            for neighbor_id, distance in neighbors:
                normalized_distance = distance / (2 * R_use)
                out.write(f"\t{neighbor_id}\t{scales.get(neighbor_id, 1.0):.2f}\t{normalized_distance:.2f}")
            out.write("\n")

    console.print(f"Saved neighbors to {output_file}", style="bold green")