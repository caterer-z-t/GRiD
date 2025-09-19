#!/usr/bin/env python3
"""
Nearest neighbors workflow for VNTR/LPA analysis

Usage:
python find_neighbors.py <batch_num> <total_batches> <max_z_range> <input_file> <output_prefix>

Example:
python find_neighbors.py 0 1 2.0 /nas/longleaf/home/catererz/work/lpa_normalized/LPA_ID_scale_zdepths.txt.gz /nas/longleaf/home/catererz/work/lpa_neighbors/LPA_find_neighbor
"""

import sys
import gzip
import time
import numpy as np
import pandas as pd
from sklearn.neighbors import NearestNeighbors

def read_compressed_file(filename: str):
    """Read gz or plain text file into pandas DataFrame"""
    opener = gzip.open if filename.endswith(".gz") else open
    return pd.read_csv(filename, sep="\t", header=None, compression="gzip" if filename.endswith(".gz") else None)

def crop_z_scores(arr: np.ndarray, z_max: float) -> np.ndarray:
    """Crop z-scores elementwise"""
    return np.clip(arr, -z_max, z_max)

def main():
    if len(sys.argv) != 6:
        print("Usage: python find_neighbors.py <batch_num> <total_batches> <max_z_range> <input_file> <output_prefix>", file=sys.stderr)
        return 1

    batch_num = int(sys.argv[1])
    total_batches = int(sys.argv[2])
    z_max = float(sys.argv[3])
    data_file = sys.argv[4]
    out_prefix = sys.argv[5]

    print(f"[INFO] Batch {batch_num}/{total_batches}, cropping to ±{z_max}")

    start_time = time.time()
    
    # --- Read input ---
    df = read_compressed_file(data_file)
    print(f"[INFO] Loaded {df.shape[0]} samples with {df.shape[1]} columns")

    # First two columns = IDs
    ids = df.iloc[:, 0].astype(str).tolist()
    scales = df.iloc[:, 1].astype(float).tolist()
    features = df.iloc[:, 2:].to_numpy(dtype=np.float32)

    # Crop z-scores
    features = crop_z_scores(features, z_max)

    N, R = features.shape
    print(f"[INFO] Data shape: {N} x {R}")

    # --- Split into batches ---
    batch_indices = [i for i in range(N) if i % total_batches == batch_num]
    X_batch = features[batch_indices, :]
    print(f"[INFO] Processing {len(batch_indices)} samples in this batch")

    # --- Fit NearestNeighbors model ---
    nn = NearestNeighbors(n_neighbors=min(501, N), metric="euclidean", n_jobs=-1)
    nn.fit(features)

    # --- Query neighbors for batch ---
    distances, indices = nn.kneighbors(X_batch)

    # --- Write output ---
    out_file = f"{out_prefix}.zMax{z_max}.txt.gz"
    with gzip.open(out_file, "wt") as fout:
        for i, row_idx in enumerate(batch_indices):
            fout.write(f"{ids[row_idx]}\t{scales[row_idx]:.2f}")
            for d, idx in zip(distances[i, 1:], indices[i, 1:]):  # skip self at [0]
                normalized_dist = d**2 / (2 * R)  # match your C++ normalization
                fout.write(f"\t{ids[idx]}\t{scales[idx]:.2f}\t{normalized_dist:.2f}")
            fout.write("\n")

    print(f"[INFO] Done in {time.time()-start_time:.2f} sec. Output: {out_file}")

if __name__ == "__main__":
    sys.exit(main())
