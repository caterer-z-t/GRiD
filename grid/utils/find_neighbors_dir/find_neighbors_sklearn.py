# grid/utils/find_neighbors_dir/find_neighbors_sklearn.py
# In[1]: Imports
import numpy as np
from sklearn.neighbors import NearestNeighbors

# In[2]: Function to find neighbors using sklearn
def find_neighbors_sklearn(data_matrix, individuals, n_neighbors=500, zmax=2.0):
    """
    Find nearest neighbors using scikit-learn.
    
    Args:
        data_matrix: np.array of shape [regions x individuals]
        individuals: list of IDs (length = number of individuals)
        n_neighbors: how many neighbors to return per individual
        zmax: max z-score clipping

    Returns:
        dict mapping individual ID to list of (neighbor ID, distance)
    """
    # Clip and replace NaNs
    X = np.clip(data_matrix.T, -zmax, zmax)  # transpose -> [individuals x regions]
    X = np.nan_to_num(X, nan=0.0)

    # Fit NearestNeighbors
    # Note: n_neighbors+1 to account for self, which we'll filter out
    nbrs = NearestNeighbors(
        n_neighbors=min(n_neighbors + 1, len(individuals)),  # Don't exceed dataset size
        algorithm="auto", 
        metric="euclidean"
    ).fit(X)

    # Get neighbors (distances and indices)
    distances, indices = nbrs.kneighbors(X)

    # Build output dictionary
    results = {}
    for i, ind in enumerate(individuals):
        # Filter out self (distance=0) and collect neighbors
        neighbor_data = [
            (individuals[j], dist) 
            for j, dist in zip(indices[i], distances[i]) 
            if j != i
        ]
        # Ensure we return exactly n_neighbors (or fewer if not enough individuals)
        results[ind] = neighbor_data[:n_neighbors]

    return results