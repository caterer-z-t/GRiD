#!/usr/bin/env python3

# In[0]: Imports
import sys
import gzip
import time
import numpy as np
import pandas as pd
import argparse
from pathlib import Path
from collections import defaultdict
from sklearn.neighbors import NearestNeighbors
import math

# Get the file path of the current script
current_file_path = Path(__file__).resolve()
current_dir = current_file_path.parent
utils_path = current_dir / "utils" / "utils.py"

# Import all contents from utils/utils.py
sys.path.insert(0, str(current_dir / "utils"))

from utils import *

# In[1]: Parse Arguments
def arg_parser():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description="Find nearest neighbors based on normalized depth data")
    parser.add_argument("--input_file", type=str, required=True, 
                       help="Input file from normalize_mosdepth (gzipped)")
    parser.add_argument("--output_prefix", type=str, required=True,
                       help="Output file prefix")
    parser.add_argument("--zmax", type=float, default=2.0,
                       help="Maximum z-score value to allow (cropping threshold)")
    parser.add_argument("--n_neighbors", type=int, default=500,
                       help="Number of nearest neighbors to output")
    parser.add_argument("--sigma2_max", type=float, default=1000.0,
                       help="Maximum variance ratio threshold")
    return parser.parse_args()

# In[1.1]: Utility Functions
def square(x):
    """Square a number"""
    return x * x

def crop(z, zmax):
    """Crop/clamp a z-score value to be within [-zmax, zmax] range"""
    return min(zmax, max(-zmax, z))

# In[2]: Read and parse the normalized data
def read_normalized_data(input_file):
    """
    Read the normalized data from your Python normalization script.
    Returns: individuals (list), regions (list), data_matrix (numpy array), scales (dict)
    
    Expected format from your normalize_mosdepth.py:
    region    ind1    ind2    ind3    ...
    start1-end1    z11    z12    z13    ...
    start2-end2    z21    z22    z23    ...
    """
    log_message(f"Reading normalized data from {input_file}")
    
    with gzip.open(input_file, 'rt') as f:
        # Read header to get individual names
        header = f.readline().strip().split('\t')
        individuals = header[1:]  # Skip 'region' column
        N = len(individuals)
        
        # Read all data rows
        regions = []
        data_rows = []
        
        for line in f:
            fields = line.strip().split('\t')
            region = fields[0]
            values = []
            
            for val_str in fields[1:]:
                if val_str == 'NA' or val_str == 'nan':
                    values.append(np.nan)
                else:
                    values.append(float(val_str))
            
            regions.append(region)
            data_rows.append(values)
    
    # Convert to numpy array: [regions x individuals]
    data_matrix = np.array(data_rows)
    R = len(regions)
    
    log_message(f"Loaded data: {N} individuals, {R} regions")
    
    # For compatibility with C++ code, create dummy scales (all 1.0)
    # In the original C++ code, scales were read from the input file
    scales = {ind: 1.0 for ind in individuals}
    
    return individuals, regions, data_matrix, scales

# In[3]: Filter regions based on variance
def filter_regions_by_variance(data_matrix, regions, sigma2_max=1000.0, frac_r=1.0):
    """
    Filter regions based on variance ratios, similar to C++ code.
    Returns: filtered indices, variance ratios
    """
    log_message("Filtering regions based on variance...")
    
    R, N = data_matrix.shape
    variance_ratios = []
    
    # Calculate variance ratio for each region
    for r in range(R):
        region_data = data_matrix[r, :]
        # Remove NaN values
        valid_data = region_data[~np.isnan(region_data)]
        
        if len(valid_data) == 0:
            variance_ratios.append(np.inf)
            continue
            
        mean_val = np.mean(valid_data)
        var_val = np.var(valid_data)
        
        # Calculate variance ratio (similar to sigma2ratio in C++)
        if mean_val != 0:
            variance_ratios.append(abs(var_val))  # Use absolute variance as proxy
        else:
            variance_ratios.append(np.inf)
    
    variance_ratios = np.array(variance_ratios)
    
    # Sort and determine thresholds (similar to C++ logic)
    sorted_vars = np.sort(variance_ratios[np.isfinite(variance_ratios)])
    if len(sorted_vars) == 0:
        log_message("Warning: No valid variance ratios found")
        return np.arange(R), variance_ratios
    
    sigma2_min = sorted_vars[int(len(sorted_vars) * (1 - frac_r))] if frac_r < 1.0 else 0
    
    # Find regions that pass filtering
    valid_indices = []
    extreme_count = 0
    
    for r in range(R):
        var_ratio = variance_ratios[r]
        if np.isfinite(var_ratio) and sigma2_min <= var_ratio <= sigma2_max:
            valid_indices.append(r)
        elif var_ratio > sigma2_max:
            extreme_count += 1
    
    valid_indices = np.array(valid_indices)
    
    log_message(f"Removed {extreme_count} of {R} regions with variance > {sigma2_max}")
    log_message(f"Keeping {len(valid_indices)} regions for analysis")
    
    return valid_indices, variance_ratios


# In[5]: Calculate distances between all individuals and batch individuals
def calculate_distances(data_matrix, valid_region_indices, batch_indices, zmax):
    """
    Calculate Euclidean distances between all individuals and batch individuals.
    Similar to C++ distance calculation logic.
    """
    log_message("Calculating distances...")
    
    R_filtered = len(valid_region_indices)
    N = data_matrix.shape[1]
    N_batch = len(batch_indices)
    
    # Extract and crop data for valid regions only
    filtered_data = data_matrix[valid_region_indices, :]  # [R_filtered x N]
    
    # Apply cropping to all data
    filtered_data = np.clip(filtered_data, -zmax, zmax)
    
    # Handle NaN values by setting them to 0 (or you could use mean imputation)
    filtered_data = np.nan_to_num(filtered_data, nan=0.0)
    
    # Extract batch data: [R_filtered x N_batch]
    batch_data = filtered_data[:, batch_indices]
    
    # Calculate squared Euclidean distances
    # distances[n, i] = sum over regions of (data[r,n] - batch_data[r,i])^2
    distances = np.zeros((N, N_batch))
    
    for i, batch_idx in enumerate(batch_indices):
        # For each individual in batch, compute distance to all individuals
        batch_profile = batch_data[:, i:i+1]  # [R_filtered x 1]
        
        # Compute squared differences: [R_filtered x N]
        squared_diffs = (filtered_data - batch_profile) ** 2
        
        # Sum over regions: [N]
        distances[:, i] = np.sum(squared_diffs, axis=0)
    
    log_message(f"Computed distances: {N} x {N_batch} matrix")
    return distances

# In[6]: Find nearest neighbors and write output
def find_neighbors_and_write(distances, individuals, scales, indices, output_prefix, zmax, n_neighbors, R_use):
    """
    Find nearest neighbors and write output in format similar to C++ version.
    """
    log_message("Finding nearest neighbors and writing output...")
    
    # Create output filename
    output_file = f"{output_prefix}.zMax{zmax:.1f}.txt.gz"
    
    N = len(individuals)
    N_batch = len(batch_individuals)
    
    with gzip.open(output_file, 'wt') as out:
        # For each individual in the batch
        for i, (batch_ind, batch_global_idx) in enumerate(zip(batch_individuals, batch_indices)):
            # Get distances for this batch individual
            dist_vector = distances[:, i]
            
            # Create (distance, individual_index) pairs
            dist_pairs = [(dist_vector[n], n) for n in range(N)]
            
            # Set self-distance to infinity to exclude from neighbors
            dist_pairs[batch_global_idx] = (float('inf'), batch_global_idx)
            
            # Sort by distance
            dist_pairs.sort(key=lambda x: x[0])
            
            # Write output: ID, scale, then top N neighbors
            out.write(f"{batch_ind}\t{scales[batch_ind]:.2f}")
            
            # Output nearest neighbors
            for j in range(min(n_neighbors, len(dist_pairs))):
                distance, neighbor_idx = dist_pairs[j]
                if distance == float('inf'):
                    continue
                    
                neighbor_id = individuals[neighbor_idx]
                neighbor_scale = scales[neighbor_id]
                normalized_distance = distance / (2 * R_use)  # Same normalization as C++
                
                out.write(f"\t{neighbor_id}\t{neighbor_scale:.2f}\t{normalized_distance:.2f}")
            
            out.write("\n")
    
    log_message(f"Output written to {output_file}")

def find_neighbors_sklearn(data_matrix, individuals, n_neighbors=500, zmax=2.0):
    """
    Find nearest neighbors using scikit-learn.
    
    Args:
        data_matrix: np.array of shape [regions x individuals]
        individuals: list of IDs (length = number of individuals)
        n_neighbors: how many neighbors to return per individual
        zmax: max z-score clipping
    """
    # Clip and replace NaNs
    X = np.clip(data_matrix.T, -zmax, zmax)  # transpose -> [individuals x regions]
    X = np.nan_to_num(X, nan=0.0)

    # Fit NearestNeighbors
    nbrs = NearestNeighbors(n_neighbors=n_neighbors+1,  # +1 to exclude self
                            algorithm="auto", 
                            metric="euclidean").fit(X)

    # Get neighbors (distances and indices)
    distances, indices = nbrs.kneighbors(X)

    # Build output dictionary
    results = {}
    for i, ind in enumerate(individuals):
        neighbor_ids = [individuals[j] for j in indices[i] if j != i]
        neighbor_dists = [dist for dist, j in zip(distances[i], indices[i]) if j != i]
        results[ind] = list(zip(neighbor_ids, neighbor_dists))

    return results

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

    log_message(f"Neighbors saved to {output_file}")


# In[7]: Main execution
def main():
    """Main function that orchestrates the neighbor finding process"""
    # Parse arguments
    args = arg_parser()
    start_time = time.time()
    
    log_message(f"Starting neighbor finding with zmax={args.zmax}")
    
    # Phase 1: Read normalized data
    individuals, regions, data_matrix, scales = read_normalized_data(args.input_file)
    
    # Phase 2: Filter regions based on variance
    valid_region_indices, variance_ratios = filter_regions_by_variance(
        data_matrix, regions, args.sigma2_max
    )
    R_use = len(valid_region_indices)
    
    # get indicies from individuals
    N = len(individuals)
    indices = list(range(N))

    # Phase 4: Calculate distances
    distances = calculate_distances(
        data_matrix, valid_region_indices, indices, args.zmax
    )
    
    # Phase 5: Find neighbors and write output
    # find_neighbors_and_write(
    #     distances, individuals, scales, indices,
    #     args.output_prefix, args.zmax, args.n_neighbors, R_use
    # )

    neighbors = find_neighbors_sklearn(
        data_matrix[valid_region_indices, :], individuals,
        n_neighbors=args.n_neighbors, zmax=args.zmax
    )

    save_neighbors(neighbors, scales, args.output_prefix, args.zmax)
    
    end_time = time.time()
    log_message(f"Neighbor finding completed in {end_time - start_time:.2f} seconds")

if __name__ == "__main__":
    main()