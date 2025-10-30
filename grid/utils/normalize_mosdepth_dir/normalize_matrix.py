# grid/utils/normalize_mosdepth_dir/normalize_matrix.py
# In[1]: Imports
import numpy as np

# In[2]: Function to Normalize Matrix
def normalize_matrix(mat):
    """
    Normalize the depth matrix in two steps:
    1) Within-individual normalization (row-wise): divide each depth by the individual's mean
    2) Across-individual normalization (column-wise): for each region, compute mean (mu) and variance (var)

    - Transform each depth x to (x - mu) / sqrt(mu) where mu > 0
    - Also compute variance ratio var/mu for each region.

    Args:
        mat: numpy array of shape (n_individuals, n_regions) with raw depths
    
    Returns:
        normalized_mat: numpy array of same shape with normalized depths
        variance_ratios: dict of {region_index: variance_ratio}
    """
    mat = mat.copy()
    # 1) within-individual (row-wise)
    row_means = np.nanmean(mat, axis=1)
    # avoid divide-by-zero
    row_means_safe = np.where(row_means == 0, np.nan, row_means)
    mat = (mat.T / row_means_safe).T  # divide each row

    # 2) across-individual
    col_means = np.nanmean(mat, axis=0)  # mu
    col_vars = np.nanvar(mat, axis=0)
    # variance ratio var/mu when mu>0
    with np.errstate(invalid="ignore", divide="ignore"):
        var_ratio = np.where(col_means > 0, col_vars / col_means, np.nan)

    # transform: (x - mu) / sqrt(mu) only where mu > 0
    mu_pos = col_means > 0
    sqrt_mu = np.sqrt(col_means, where=mu_pos, out=np.full_like(col_means, np.nan))
    # broadcast: subtract mu then divide
    # for columns with mu_pos False, leave as is (or NaN)
    mat[:, mu_pos] = (mat[:, mu_pos] - col_means[mu_pos]) / sqrt_mu[mu_pos]

    # convert var_ratio to python dict of index->value for ease of selection later
    variance_ratios = {i: float(var_ratio[i]) for i in range(len(col_means)) if not np.isnan(var_ratio[i])}

    return mat, variance_ratios
