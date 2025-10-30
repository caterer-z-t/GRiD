# grid/utils/normalize_mosdepth_dir/write_normalized_output.py
# In[1]: Imports
import gzip
import numpy as np
from pathlib import Path

# In[2]: Function to Write Normalized Output
def write_normalized_output(
    mat: np.ndarray,
    individuals_order: list[str],
    regions_list: list[tuple[int, int]],
    selected_indices: list[int],
    output_file: Path,
):
    """
    Write normalized matrix to gzipped TSV.

    Args:
        mat: normalized matrix (n_individuals × n_regions)
        individuals_order: list of sample IDs
        regions_list: list of region tuples
        selected_indices: list of region indices to keep
        output_file: Path to write .tsv.gz file
    """
    with gzip.open(output_file, "wt") as out:
        out.write("region\t" + "\t".join(individuals_order) + "\n")
        for j in selected_indices:
            start, end = regions_list[j]
            vals = [
                "NA" if np.isnan(mat[i, j]) else f"{mat[i, j]:.6g}"
                for i in range(len(individuals_order))
            ]
            out.write(f"{start}-{end}\t" + "\t".join(vals) + "\n")
