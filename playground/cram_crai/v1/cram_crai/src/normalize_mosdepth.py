#!/usr/bin/env python3
# In[1]: Imports
import sys
import gzip
import math
import numpy as np
import time
from collections import defaultdict
import os
import argparse
from pathlib import Path

# Get the file path of the current script
current_file_path = Path(__file__).resolve()
current_dir = current_file_path.parent
utils_path = current_dir / "utils" / "utils.py"

# Import all contents from utils/utils.py
sys.path.insert(0, str(current_dir / "utils"))

from utils import *

# In[2]: 1. Preparation
def arg_parser():
    parser = argparse.ArgumentParser(description="Normalize mosdepth data")
    parser.add_argument("--mosdepth_dir", type=str, required=True, help="Directory containing mosdepth files")
    parser.add_argument("--repeat_mask", type=str, required=True, help="Path to repeat mask BED file")
    parser.add_argument("--chromosome", type=str, required=True, help="Chromosome to analyze")
    parser.add_argument("--start", type=int, required=True, help="Start position")
    parser.add_argument("--end", type=int, required=True, help="End position")
    parser.add_argument("--threads", type=int, default=1, help="Number of threads to use")
    parser.add_argument("--min_depth", type=int, default=20, help="Minimum depth threshold")
    parser.add_argument("--max_depth", type=int, default=100, help="Maximum depth threshold")
    parser.add_argument("--output_file", type=str, default="normalized_output.txt.gz", help="Output file path")
    return parser.parse_args()


# In[3]: Decide which regions to extract
def get_individuals(mosdepth_dir: str) -> dict[str, Path]:
    """Map individual IDs to their mosdepth.global.dist file paths."""
    mosdepth_dir = Path(mosdepth_dir)
    return {
        f.stem.split(".")[0]: f
        for f in mosdepth_dir.glob("*.mosdepth.global.dist.txt")
    }


def regions_bed_gz(global_dist_file: Path) -> Path:
    """Derive the .regions.bed.gz file path from a global.dist file path."""
    return global_dist_file.with_name(
        global_dist_file.name.replace(".mosdepth.global.dist.txt", ".regions.bed.gz")
    )


def extract_reasonable_depth_regions(
    individuals: dict[str, Path],
    chromosome: str,
    start: int,
    end: int,
    min_depth: int = 20,
    max_depth: int = 100,
) -> dict[str, list[tuple[int, int]]]:
    """Extract regions within depth bounds for each individual."""
    regions_to_extract: dict[str, list[tuple[int, int, float]]] = defaultdict(list)

    for individual_id, global_dist_file in individuals.items():
        bed_gz = regions_bed_gz(global_dist_file)

        log_message(f"Processing {bed_gz} for individual {individual_id}")
        
        if not bed_gz.exists():
            log_message(f"Warning: {bed_gz} does not exist.")
            continue

        with gzip.open(bed_gz, "rt") as f:
            for line in f:
                if not line.startswith(f"chr{chromosome}"):
                    continue

                fields = line.strip().split("\t")
                if len(fields) < 4:
                    continue  # skip malformed lines

                chrom, region_start, region_end, depth = fields[0], int(fields[1]), int(fields[2]), float(fields[3])

                # Only skip if NEITHER condition is true
                if depth <= 0 and (region_end < start or region_start > end):
                    continue

                if depth < min_depth or depth > max_depth:
                    continue
                
                # Store if you want to return later
                regions_to_extract[individual_id].append((region_start, region_end, depth))

    return regions_to_extract

def extract_reasonable_depth_regions_excluded(
    individuals: dict[str, Path],
    chromosome: str,
    start: int,
    end: int,
    min_depth: int = 20,
    max_depth: int = 100,
    excluded: dict[str, set[int]] = None
) -> dict[str, list[tuple[int, int, float]]]:
    """Extract regions within depth bounds for each individual, excluding repeats."""
    regions_to_extract: dict[str, list[tuple[int, int, float]]] = defaultdict(list)

    for individual_id, global_dist_file in individuals.items():
        bed_gz = regions_bed_gz(global_dist_file)

        log_message(f"Processing {bed_gz} for individual {individual_id}")
        
        if not bed_gz.exists():
            log_message(f"Warning: {bed_gz} does not exist.")
            continue

        with gzip.open(bed_gz, "rt") as f:
            for line in f:
                if not line.startswith(f"chr{chromosome}"):
                    continue

                fields = line.strip().split("\t")
                if len(fields) < 4:
                    continue  # skip malformed lines

                chrom, region_start, region_end, depth = fields[0], int(fields[1]), int(fields[2]), float(fields[3])

                # Condition 1 OR 2: depth > 0 OR overlap with [start, end]
                if not (depth > 0 or (region_end >= start and region_start <= end)):
                    continue

                # Repeat filtering
                if excluded is not None:
                    region_kb = set(range(region_start // 1000, region_end // 1000 + 1))
                    if region_kb & excluded.get(chrom, set()):
                        # log_message(f"Skipping repeat-overlapping region {chrom}:{region_start}-{region_end}")
                        continue

                # Passed filters → log the line
                regions_to_extract[individual_id].append((region_start, region_end, depth))
                
    return regions_to_extract

# In[4]: Exclude regions overlapping repeats/VNTR

def load_repeat_mask(repeat_bed: str) -> dict[str, set[int]]:
    """
    Read repeat/VNTR regions into a dictionary of excluded kb bins per chromosome.
    Returns: {chrom: set(kb_indices)}
    """
    excluded = defaultdict(set)
    with open(repeat_bed) as f:
        for line in f:
            if line.startswith("#"):  # skip header/comment
                continue
            chrom, start, end = line.strip().split()[:3]
            start, end = int(start), int(end)
            for kb in range(start // 1000, end // 1000 + 1):
                excluded[chrom].add(kb)
    return excluded

def build_depth_matrix(regions_to_extract: dict[str, list[tuple[int, int, float]]]):
    """
    Build depth matrix: {individual: {region: depth}} 
    directly from regions_to_extract, which already has (start, end, depth).
    """
    depth_matrix = defaultdict(dict)

    for individual_id, regions in regions_to_extract.items():
        for start, end, depth in regions:
            depth_matrix[individual_id][(start, end)] = depth

    return depth_matrix

def normalize_within_individual(depth_matrix):
    """
    Normalize depths within each individual by their mean depth.
    """
    for ind, regions in depth_matrix.items():
        if not regions:
            continue
        mean_depth = np.mean(list(regions.values()))
        for r in regions:
            regions[r] /= mean_depth
    return depth_matrix

def normalize_across_individuals(depth_matrix):
    """
    Normalize per-region across all individuals: z = (x - mu) / sqrt(mu).
    Returns normalized matrix and variance ratios.
    """
    # collect all regions
    all_regions = set()
    for regions in depth_matrix.values():
        all_regions.update(regions.keys())

    variance_ratios = {}
    for region in all_regions:
        values = [depth_matrix[ind].get(region, np.nan) for ind in depth_matrix]
        values = np.array([v for v in values if not np.isnan(v)])
        if len(values) == 0:
            continue

        mu = values.mean()
        var = values.var()
        if mu > 0:
            var_ratio = var / mu
            variance_ratios[region] = var_ratio

            # transform each individual's value
            for ind in depth_matrix:
                if region in depth_matrix[ind]:
                    depth_matrix[ind][region] = (depth_matrix[ind][region] - mu) / math.sqrt(mu)

    return depth_matrix, variance_ratios

def select_high_variance_regions(variance_ratios, top_frac=0.1):
    """
    Select top fraction of regions by variance/mean ratio.
    """
    n_keep = max(1, int(len(variance_ratios) * top_frac))
    sorted_regions = sorted(variance_ratios.items(), key=lambda x: x[1], reverse=True)
    return dict(sorted_regions[:n_keep])

def write_output(depth_matrix, variance_ratios, output_file):
    individuals = list(depth_matrix.keys())
    regions = sorted(variance_ratios.keys())

    with gzip.open(output_file, "wt") as out:
        # header
        out.write("region\t" + "\t".join(individuals) + "\n")

        # per-region normalized depths
        for region in regions:
            line = [f"{region[0]}-{region[1]}"]
            for ind in individuals:
                val = depth_matrix[ind].get(region, "NA")
                line.append(str(val))
            out.write("\t".join(line) + "\n")


# In[]: Main execution
if __name__ == "__main__":
    # 1. Preparation
    args = arg_parser()
    start_time = time.time()

    # Get individuals and their mosdepth files
    individuals = get_individuals(args.mosdepth_dir)
    log_message(f"Found {len(individuals)} individuals.")

    # run on first 5 individuals for testing
    # individuals = dict(list(individuals.items())[:5])
    # log_message(f"Processing {len(individuals)} individuals for testing.")

    # Load repeat mask
    excluded = load_repeat_mask(args.repeat_mask)

    # Phase 1 + 2: extract filtered regions
    regions_to_extract = extract_reasonable_depth_regions_excluded(
        individuals,
        args.chromosome,
        args.start,
        args.end,
        min_depth=args.min_depth,
        max_depth=args.max_depth,
        excluded=excluded
    )

    # Phase 3: build depth matrix
    depth_matrix = build_depth_matrix(regions_to_extract)
    
    # Phase 4: within-individual normalization
    depth_matrix = normalize_within_individual(depth_matrix)

    # Phase 5: across-individual normalization
    depth_matrix, variance_ratios = normalize_across_individuals(depth_matrix)

    # Phase 6: keep high variance regions
    variance_ratios = select_high_variance_regions(variance_ratios, top_frac=0.1)

    # Phase 7: write results
    write_output(depth_matrix, variance_ratios, args.output_file)

    end_time = time.time()
    log_message(f"Normalization completed in {end_time - start_time:.2f} seconds.")