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
from rich.console import Console

console = Console()

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
        
        if not bed_gz.exists():
            console.log(f"[yellow]Warning: {bed_gz} does not exist for {individual_id}[/yellow]")
            continue

        region_count = 0
        with gzip.open(bed_gz, "rt") as f:
            for line in f:
                # Match chromosome (with or without 'chr' prefix)
                if chromosome.startswith('chr'):
                    if not line.startswith(chromosome):
                        continue
                else:
                    if not line.startswith(f"chr{chromosome}") and not line.startswith(chromosome):
                        continue

                fields = line.strip().split("\t")
                if len(fields) < 4:
                    continue  # skip malformed lines

                chrom = fields[0]
                region_start = int(fields[1])
                region_end = int(fields[2])
                depth = float(fields[3])

                # Check if region overlaps with our target range [start, end]
                if region_end < start or region_start > end:
                    continue

                # Apply depth filters
                if depth < min_depth or depth > max_depth:
                    continue

                # Repeat filtering - check if region overlaps with excluded kb bins
                if excluded is not None:
                    region_kb_bins = set(range(region_start // 1000, region_end // 1000 + 1))
                    if region_kb_bins & excluded.get(chrom, set()):
                        continue

                # Passed all filters
                regions_to_extract[individual_id].append((region_start, region_end, depth))
                region_count += 1
        
        if region_count > 0:
            console.log(f"[green]Found {region_count} regions for {individual_id}[/green]")
        else:
            console.log(f"[yellow]No regions found for {individual_id}[/yellow]")
                
    return regions_to_extract

# In[4]: Exclude regions overlapping repeats/VNTR

def load_repeat_mask(repeat_bed: str, chromosome: str = None) -> dict[str, set[int]]:
    """
    Read repeat/VNTR regions into a dictionary of excluded kb bins per chromosome.
    Returns: {chrom: set(kb_indices)}
    """
    excluded = defaultdict(set)
    console.log(f"Loading repeat mask from {repeat_bed}")
    
    with open(repeat_bed) as f:
        for line in f:
            if line.startswith("#"):  # skip header/comment
                continue
            parts = line.strip().split()
            if len(parts) < 3:
                continue
                
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            
            # Filter by chromosome if specified
            if chromosome is not None:
                target_chr = chromosome if chromosome.startswith('chr') else f'chr{chromosome}'
                if chrom != target_chr:
                    continue
            
            # Mark kb bins that overlap this repeat region
            for kb in range(start // 1000, end // 1000 + 1):
                excluded[chrom].add(kb)
    
    total_excluded = sum(len(bins) for bins in excluded.values())
    console.log(f"Loaded {total_excluded} excluded kb bins across {len(excluded)} chromosomes")
    
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
        if mean_depth > 0:
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

    console.log(f"Normalizing across {len(all_regions)} unique regions")

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
    console.log(f"Selecting top {n_keep} regions (top {top_frac*100}%) by variance")
    return dict(sorted_regions[:n_keep])

def write_output(depth_matrix, variance_ratios, output_file):
    individuals = list(depth_matrix.keys())
    regions = sorted(variance_ratios.keys())

    console.log(f"Writing {len(regions)} regions x {len(individuals)} individuals to {output_file}")

    with gzip.open(output_file, "wt") as out:
        # header
        out.write("region\t" + "\t".join(individuals) + "\n")

        # per-region normalized depths
        for region in regions:
            line = [f"{region[0]}-{region[1]}"]
            for ind in individuals:
                val = depth_matrix[ind].get(region, "NA")
                if val != "NA":
                    line.append(f"{val:.3f}")
                else:
                    line.append("NA")
            out.write("\t".join(line) + "\n")

def run_normalize_mosdepth(
    mosdepth_dir,
    output_file,
    repeat_mask,
    chrom,
    start,
    end,
    min_depth=20,
    max_depth=100,
    top_frac=0.1,
    threads=1,
):

    # 1. Preparation
    start_time = time.time()

    console.log(f"[bold]Starting normalization pipeline[/bold]")
    console.log(f"Target region: chr{chrom}:{start}-{end}")
    console.log(f"Depth filters: {min_depth}-{max_depth}")

    # Get individuals and their mosdepth files
    individuals = get_individuals(mosdepth_dir)
    console.log(f"Found {len(individuals)} individuals")

    if len(individuals) == 0:
        console.log("[red]ERROR: No individuals found! Check mosdepth_dir path.[/red]")
        return

    # Load repeat mask
    excluded = load_repeat_mask(repeat_mask, chromosome=chrom)

    # Phase 1 + 2: extract filtered regions
    regions_to_extract = extract_reasonable_depth_regions_excluded(
        individuals,
        chrom,
        start,
        end,
        min_depth=min_depth,
        max_depth=max_depth,
        excluded=excluded
    )

    console.log(f"Extracted regions for {len(regions_to_extract)} individuals")
    
    if len(regions_to_extract) == 0:
        console.log("[red]ERROR: No regions extracted! Check your filters and data.[/red]")
        return

    # Phase 3: build depth matrix
    depth_matrix = build_depth_matrix(regions_to_extract)
    
    total_regions = sum(len(regions) for regions in depth_matrix.values())
    console.log(f'Built depth matrix: {len(depth_matrix)} individuals, {total_regions} total region entries')

    # Phase 4: within-individual normalization
    depth_matrix = normalize_within_individual(depth_matrix)
    console.log("Completed within-individual normalization")

    # Phase 5: across-individual normalization
    depth_matrix, variance_ratios = normalize_across_individuals(depth_matrix)
    console.log(f"Completed across-individual normalization: {len(variance_ratios)} regions with variance")

    # Phase 6: keep high variance regions
    variance_ratios = select_high_variance_regions(variance_ratios, top_frac=top_frac)

    # Phase 7: write results
    write_output(depth_matrix, variance_ratios, output_file)

    end_time = time.time()
    console.log(f"[bold green]Normalization completed in {end_time - start_time:.2f} seconds[/bold green]")

# In[]: Main execution
if __name__ == "__main__":

    args = arg_parser()

    run_normalize_mosdepth(
        mosdepth_dir=args.mosdepth_dir,
        output_file=args.output_file,
        repeat_mask=args.repeat_mask,
        chrom=args.chromosome,
        start=args.start,
        end=args.end,
        min_depth=args.min_depth,
        max_depth=args.max_depth,
        top_frac=0.1,
        threads=args.threads
    )