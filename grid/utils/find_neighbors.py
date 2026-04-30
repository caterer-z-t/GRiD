# In[0]: Imports
from pathlib import Path
import numpy as np
from sklearn.neighbors import NearestNeighbors
import gzip

from .utils import log, progress_bar


# In[1]: Main function to find neighbors
def find_neighbors(config, console):
    """
    Find nearest neighbors from normalized coverage data for all samples.

    Mirrors the C++ find_neighbors.cpp logic:
      1. Read N, R, per-region means (header row 0) — means unused here but
         read to advance the file pointer.
      2. Read per-region sigma2 ratios (header row 1) — used to filter regions.
      3. Clip z-scores to [-zmax, zmax] and replace NaN → 0.
      4. Filter regions: keep those with sigma2_min <= ratio <= sigma2_max.
      5. Compute pairwise Euclidean distances on the filtered z-score matrix.
      6. For each individual output the top N_OUTPUT nearest neighbors (excluding
         self), with normalised distance = raw_dist / (2 * R_use).

    Args:
        config:  configuration dictionary
        console: Rich console for logging
    """
    try:
        threads = config.get("threads", 1)
        zmax = config["mosdepth"]["neighbors"].get("zmax", 2.0)
        sigma2_max = config["mosdepth"]["neighbors"].get("sigma2_max", 1000.0)
        n_neighbors = config["mosdepth"]["neighbors"].get("num_neighbors", 500)  # C++ outputs 500
        frac_r = config["mosdepth"]["neighbors"].get("frac_r", 1.0)  # C++ hardcodes 1

        input_file_prefix = config["mosdepth"]["normalize"].get("output_file_prefix")
        output_file_type = config.get("output_file_type", "tsv")
        output_dir = config.get("output_dir", ".")
        input_file = Path(f"{output_dir}/{input_file_prefix}.{output_file_type}.gz")

        output_prefix = config["mosdepth"]["neighbors"].get(
            "output_file_prefix", "neighbor_coverage"
        )
        # C++ writes one file per batch; Python processes all at once → single output file
        output_file = Path(output_dir) / f"{output_prefix}.zMax{zmax:.1f}.{output_file_type}.gz"
    except Exception as e:
        log(console, f"Config error: {e}", style="danger")
        return

    output_file.parent.mkdir(parents=True, exist_ok=True)

    # --- Step 1 & 2: read file ---
    individuals, sigma2ratios, data_matrix, scales = read_normalized_data(input_file)
    N, R = data_matrix.shape  # [individuals x regions]

    # --- Step 3: clip z-scores and replace NaN → 0  ---
    clipped = np.clip(data_matrix, -zmax, zmax)
    clipped = np.nan_to_num(clipped, nan=0.0)

    # --- Step 4: filter regions by sigma2 ratio ---
    valid_indices, R_use = filter_regions_by_variance(
        sigma2ratios, frac_r=frac_r, sigma2_max=sigma2_max, console=console
    )

    filtered = clipped[:, valid_indices]  # [individuals x R_use]

    # --- Step 5 & 6: find neighbors and write output ---
    with progress_bar(console, total=N, description="Finding neighbors...") as (progress, task):
        neighbors = find_neighbors_sklearn(
            filtered,
            individuals,
            n_neighbors=n_neighbors,
        )
        progress.advance(task, N)

    save_neighbors(neighbors, scales, output_file, zmax, R_use)
    log(console, f"Saved neighbors to {output_file}", style="success")


# In[2]: Read normalized data
def read_normalized_data(input_file: Path):
    """
    Parse the gzipped normalized depth file produced by normalize_mosdepth.

    File layout (matches write_normalized_output):
      Line 0 : N <tab> Rwant <tab> mu_1 ... mu_Rwant
      Line 1 : N <tab> Rwant <tab> varRatio_1 ... varRatio_Rwant
      Line 2+: ID <tab> scale <tab> z_1 ... z_Rwant

    Args:
        input_file: Path to the .tsv.gz output of normalize_mosdepth

    Returns:
        individuals  : list of sample IDs (length N)
        sigma2ratios : np.ndarray of per-region variance ratios (length Rwant)
        data_matrix  : np.ndarray shape [N x Rwant] of z-scores
        scales       : dict {sample_id: scale_value}
    """
    individuals = []
    scales = {}
    rows = []
    sigma2ratios = None

    with gzip.open(input_file, "rt") as f:
        # Header row 0: N, Rwant, mu_1 ... mu_Rwant  (means — read but not used here)
        _ = f.readline()

        # Header row 1: N, Rwant, varRatio_1 ... varRatio_Rwant
        parts = f.readline().strip().split("\t")
        # parts[0]=N, parts[1]=Rwant, parts[2:]=ratios
        sigma2ratios = np.array([np.nan if v in ("NA", "nan") else float(v) for v in parts[2:]])

        # Data rows: ID, scale, z_1 ... z_Rwant
        for line in f:
            parts = line.strip().split("\t")
            ind_id = parts[0]
            scale = float(parts[1])
            zvals = [np.nan if v in ("NA", "nan") else float(v) for v in parts[2:]]
            individuals.append(ind_id)
            scales[ind_id] = scale
            rows.append(zvals)

    data_matrix = np.array(rows, dtype=float)  # [N x Rwant]
    return individuals, sigma2ratios, data_matrix, scales


# In[3]: Filter regions by variance ratio
def filter_regions_by_variance(
    sigma2ratios: np.ndarray,
    frac_r: float = 1.0,
    sigma2_max: float = 1000.0,
    console=None,
) -> tuple[np.ndarray, int]:
    """
    Select region indices whose sigma2 ratio passes the [sigma2min, sigma2max]
    filter, exactly mirroring the C++ logic.

    Args:
        sigma2ratios : per-region variance ratios from the file header (length R)
        frac_r       : fraction of regions to retain by lower bound (C++ = 1.0)
        sigma2_max   : upper bound on variance ratio (C++ = 1000)
        console      : Rich console for logging

    Returns:
        valid_indices : np.ndarray of int indices that pass the filter
        R_use         : number of retained regions
    """
    R = len(sigma2ratios)

    # Use only finite ratios to compute the lower-bound threshold
    finite_mask = np.isfinite(sigma2ratios)
    finite_vals = np.sort(sigma2ratios[finite_mask])

    if len(finite_vals) == 0:
        log(console, "Warning: no finite variance ratios — keeping all regions", style="warning")
        return np.arange(R), R

    # When frac_r=1.0 this index is 0, so sigma2min = the smallest finite ratio
    lower_idx = int(R * (1.0 - frac_r))
    lower_idx = min(lower_idx, len(finite_vals) - 1)
    sigma2_min = float(finite_vals[lower_idx])

    extreme_count = int(np.sum(sigma2ratios > sigma2_max))
    if extreme_count and console:
        log(
            console,
            f"Removed {extreme_count} / {R} regions with sigma2ratio > {sigma2_max}",
            style="warning",
        )

    valid_mask = finite_mask & (sigma2ratios >= sigma2_min) & (sigma2ratios <= sigma2_max)
    valid_indices = np.where(valid_mask)[0]
    R_use = len(valid_indices)

    return valid_indices, R_use


# In[4]: Find neighbors using sklearn
def find_neighbors_sklearn(
    data_matrix: np.ndarray,
    individuals: list[str],
    n_neighbors: int = 500,
) -> dict[str, list[tuple[str, float]]]:
    """
    Find nearest neighbors using scikit-learn's NearestNeighbors.

    The C++ computes raw squared Euclidean distances then writes
    dist / (2 * R_use).  sklearn returns Euclidean (not squared), so we
    square the distances before normalising to match — the normalisation
    itself happens in save_neighbors.

    Clipping and NaN-filling must be done BEFORE calling this function
    (handled in find_neighbors).

    Args:
        data_matrix : np.ndarray [N x R_use], already clipped and NaN-filled
        individuals : list of N sample IDs
        n_neighbors : number of neighbors to return per individual (C++ = 500)

    Returns:
        dict {individual_id: [(neighbor_id, squared_euclidean_distance), ...]}
        Length of each list = min(n_neighbors, N-1).
    """
    N = len(individuals)
    k = min(n_neighbors + 1, N)  # +1 because sklearn includes self

    nbrs = NearestNeighbors(
        n_neighbors=k,
        algorithm="auto",
        metric="euclidean",
    ).fit(data_matrix)

    distances, indices = nbrs.kneighbors(data_matrix)

    results = {}
    for i, ind in enumerate(individuals):
        neighbor_list = []
        for j, dist in zip(indices[i], distances[i]):
            if j == i:
                continue  # exclude self (C++ sets dist[n_i] = 1e9 then sorts)
            # Store squared distance to match C++ accumulation of sq(z - zs[r])
            neighbor_list.append((individuals[j], dist**2))
            if len(neighbor_list) == n_neighbors:
                break
        results[ind] = neighbor_list

    return results


# In[5]: Save neighbors
def save_neighbors(
    neighbors_dict: dict[str, list[tuple[str, float]]],
    scales: dict[str, float],
    output_file: Path,
    zmax: float,
    R_use: int,
) -> None:
    """
    Write neighbors to a gzipped text file.

    Output format per line (matches C++):
        ID  scale  nb1_ID  nb1_scale  nb1_norm_dist  nb2_ID  ...

    where norm_dist = squared_euclidean_distance / (2 * R_use),
    exactly as in C++: fout << dists[n*Nbatch+i] / (2*Ruse)

    NOTE: the original Python constructed the output filename here again,
    causing a double-suffix bug (.zMax2.0.txt.gz.zMax2.0.txt.gz).
    The caller now passes the fully-resolved Path directly.

    Args:
        neighbors_dict : {id: [(neighbor_id, squared_dist), ...]}
        scales         : {id: scale_value}
        output_file    : fully-resolved output Path (no filename construction here)
        zmax           : z-score clip value (unused in output, kept for logging)
        R_use          : number of regions used, for distance normalisation
    """
    if R_use == 0:
        R_use = 1  # guard against division by zero

    with gzip.open(output_file, "wt") as out:
        for ind, neighbors in neighbors_dict.items():
            line = f"{ind}\t{scales.get(ind, 1.0):.2f}"
            for neighbor_id, sq_dist in neighbors:
                norm_dist = sq_dist / (2 * R_use)
                line += f"\t{neighbor_id}\t{scales.get(neighbor_id, 1.0):.2f}\t{norm_dist:.2f}"
            out.write(line + "\n")
