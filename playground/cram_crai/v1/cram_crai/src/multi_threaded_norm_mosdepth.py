#!/usr/bin/env python3
# optimized_parallel_mosdepth.py
import sys
import gzip
import math
import numpy as np
import time
from collections import defaultdict
import os
import argparse
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

# Insert utils path (same as your original)
current_file_path = Path(__file__).resolve()
current_dir = current_file_path.parent
sys.path.insert(0, str(current_dir / "utils"))
from utils import *  # noqa: E402,F401

# ---------------------------
# Argument parser (unchanged)
# ---------------------------
def arg_parser():
    parser = argparse.ArgumentParser(description="Normalize mosdepth data")
    parser.add_argument("--mosdepth_dir", type=str, required=True, help="Directory containing mosdepth files")
    parser.add_argument("--repeat_mask", type=str, required=True, help="Path to repeat mask BED file")
    parser.add_argument("--chromosome", type=str, required=True, help="Chromosome to analyze")
    parser.add_argument("--start", type=int, required=True, help="Start position")
    parser.add_argument("--end", type=int, required=True, help="End position")
    parser.add_argument("--threads", type=int, default=1, help="Number of threads/processes to use")
    parser.add_argument("--min_depth", type=int, default=20, help="Minimum depth threshold")
    parser.add_argument("--max_depth", type=int, default=100, help="Maximum depth threshold")
    parser.add_argument("--output_file", type=str, default="normalized_output.txt.gz", help="Output file path")
    return parser.parse_args()

# ---------------------------
# Helpers
# ---------------------------
def get_individuals(mosdepth_dir: str) -> dict[str, Path]:
    mosdepth_dir = Path(mosdepth_dir)
    return {
        f.stem.split(".")[0]: f
        for f in mosdepth_dir.glob("*.mosdepth.global.dist.txt")
    }

def regions_bed_gz(global_dist_file: Path) -> Path:
    return global_dist_file.with_name(
        global_dist_file.name.replace(".mosdepth.global.dist.txt", ".regions.bed.gz")
    )

# Worker function run in parallel per individual
def _process_one_individual(args_tuple):
    """
    This runs in its own process. Returns (individual_id, [(start,end,depth),...])
    args_tuple contains:
      individual_id, global_dist_file, chromosome, start, end, min_depth, max_depth, excluded
    """
    (individual_id, global_dist_file, chromosome, start, end, min_depth, max_depth, excluded) = args_tuple
    bed_gz = regions_bed_gz(global_dist_file)
    results = []
    if not bed_gz.exists():
        # We cannot call your main process log_message from here reliably — return empty and main will log
        return individual_id, results

    try:
        with gzip.open(bed_gz, "rt") as f:
            for line in f:
                if not line.startswith(f"chr{chromosome}"):
                    continue
                fields = line.strip().split("\t")
                if len(fields) < 4:
                    continue
                chrom = fields[0]
                region_start = int(fields[1])
                region_end = int(fields[2])
                depth = float(fields[3])

                # Keep if depth>0 or overlaps target window
                if not (depth > 0 or (region_end >= start and region_start <= end)):
                    continue

                # depth filter
                if depth < min_depth or depth > max_depth:
                    continue

                # repeat exclusion (kb bins)
                if excluded is not None:
                    region_kb = set(range(region_start // 1000, region_end // 1000 + 1))
                    if region_kb & excluded.get(chrom, set()):
                        continue

                results.append((region_start, region_end, depth))
    except Exception:
        # If file corrupt or other IO error, just return empty (main can log)
        return individual_id, []

    return individual_id, results

# Load repeat mask (unchanged but faster loop)
def load_repeat_mask(repeat_bed: str) -> dict[str, set[int]]:
    excluded = defaultdict(set)
    with open(repeat_bed) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split()
            if len(parts) < 3:
                continue
            chrom, start_str, end_str = parts[0], parts[1], parts[2]
            try:
                start, end = int(start_str), int(end_str)
            except ValueError:
                continue
            for kb in range(start // 1000, end // 1000 + 1):
                excluded[chrom].add(kb)
    return excluded

# ---------------------------
# Vectorized building & normalization
# ---------------------------
def build_matrix_from_regions(regions_to_extract, individuals_order=None):
    """
    Build:
      - regions_list: list of region tuples (start,end)
      - inds: list of individual ids
      - mat: numpy array shape (n_inds, n_regions) with np.nan for missing
    """
    if individuals_order is None:
        individuals_order = sorted(regions_to_extract.keys())
    # Collect all regions as tuples (start,end)
    region_set = set()
    for ind in individuals_order:
        for s, e, d in regions_to_extract.get(ind, []):
            region_set.add((s, e))
    regions_list = sorted(region_set)

    n_inds = len(individuals_order)
    n_regs = len(regions_list)
    mat = np.full((n_inds, n_regs), np.nan, dtype=float)

    region_index = {r: i for i, r in enumerate(regions_list)}
    ind_index = {ind: i for i, ind in enumerate(individuals_order)}

    # Fill matrix with raw depths
    for ind, regions in regions_to_extract.items():
        i = ind_index[ind]
        for s, e, d in regions:
            j = region_index.get((s, e))
            if j is not None:
                mat[i, j] = d

    return individuals_order, regions_list, mat

def normalize_matrix(mat):
    """
    1) Within-individual normalization: divide each row by the row mean (ignoring nan).
    2) Across-individual normalization: for each column (region),
       compute mu = mean(column ignoring nan), var = var(column ignoring nan),
       variance_ratio = var/mu if mu > 0 else np.nan
       transform entries x -> (x - mu) / sqrt(mu) (only if mu>0); NaNs remain.
    Returns transformed matrix and variance_ratios dict {region_index: ratio}.
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

def select_top_regions(variance_ratios: dict, regions_list, top_frac=0.1):
    """
    variance_ratios: {region_index: ratio}
    returns set of selected region indices
    """
    if not variance_ratios:
        return set()
    items = sorted(variance_ratios.items(), key=lambda x: x[1], reverse=True)
    n_keep = max(1, int(len(items) * top_frac))
    selected = set(idx for idx, _ in items[:n_keep])
    # map back to region tuples
    selected_regions = [regions_list[i] for i in sorted(selected)]
    return selected, selected_regions

def write_output_from_matrix(output_file, individuals_order, regions_list, mat, selected_region_indices):
    """
    Write gzipped table: header region\tind1\tind2...
    Only write selected regions in order.
    """
    with gzip.open(output_file, "wt") as out:
        out.write("region\t" + "\t".join(individuals_order) + "\n")
        for j in sorted(selected_region_indices):
            region = regions_list[j]
            row_vals = []
            for i in range(len(individuals_order)):
                val = mat[i, j]
                if np.isnan(val):
                    row_vals.append("NA")
                else:
                    # Format float with reasonable precision
                    row_vals.append(f"{val:.6g}")
            out.write(f"{region[0]}-{region[1]}\t" + "\t".join(row_vals) + "\n")

# ---------------------------
# Main
# ---------------------------
if __name__ == "__main__":
    args = arg_parser()
    start_time = time.time()

    # Get individuals
    individuals = get_individuals(args.mosdepth_dir)
    log_message(f"Found {len(individuals)} individuals.")

    if len(individuals) == 0:
        log_message("No individuals found. Exiting.")
        sys.exit(1)

    # Load repeat mask
    excluded = load_repeat_mask(args.repeat_mask)
    log_message("Loaded repeat mask.")

    # Prepare tasks for parallel processing
    tasks = []
    for ind_id, global_file in individuals.items():
        tasks.append((ind_id, global_file, args.chromosome, args.start, args.end, args.min_depth, args.max_depth, excluded))

    # Run parallel extraction
    regions_to_extract = defaultdict(list)
    missing_files = []
    n_workers = max(1, args.threads)

    log_message(f"Starting parallel extraction with {n_workers} workers...")
    with ProcessPoolExecutor(max_workers=n_workers) as exe:
        futures = {exe.submit(_process_one_individual, t): t[0] for t in tasks}
        for fut in as_completed(futures):
            ind = futures[fut]
            try:
                ind_id, result = fut.result()
                if not result:
                    # file missing or empty result; log if missing file exists?
                    # check file path quick
                    global_file = individuals[ind_id]
                    if not regions_bed_gz(global_file).exists():
                        log_message(f"Warning: regions file for {ind_id} missing: {regions_bed_gz(global_file)}")
                    # else it's just empty after filtering
                else:
                    regions_to_extract[ind_id] = result
            except Exception as e:
                log_message(f"Error processing {ind}: {e}")

    n_regions_found = sum(len(v) for v in regions_to_extract.values())
    log_message(f"Extracted {n_regions_found} region entries across {len(regions_to_extract)} individuals.")

    # Build matrix
    individuals_order = sorted(regions_to_extract.keys())
    if not individuals_order:
        log_message("No regions to process after extraction. Exiting.")
        sys.exit(1)

    individuals_order, regions_list, mat = build_matrix_from_regions(regions_to_extract, individuals_order)
    log_message(f"Built matrix: {mat.shape[0]} individuals x {mat.shape[1]} regions.")

    # Normalize
    mat_normed, variance_ratios = normalize_matrix(mat)
    log_message("Completed normalization (within-individual and across-individual).")

    # Select top variance regions
    selected_indices, selected_regions = select_top_regions(variance_ratios, regions_list, top_frac=0.1)
    log_message(f"Selected {len(selected_indices)} high-variance regions (top 10%).")

    # Write output
    write_output_from_matrix(args.output_file, individuals_order, regions_list, mat_normed, selected_indices)
    end_time = time.time()
    log_message(f"Normalization completed in {end_time - start_time:.2f} seconds.")
