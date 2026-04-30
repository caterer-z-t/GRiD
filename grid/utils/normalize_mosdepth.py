# In[0]: Imports
from pathlib import Path
import glob
import gzip
from collections import defaultdict
import numpy as np

from .utils import (
    log,
    get_samples,
    # create_region_string,
    setup_output_file,
    progress_bar,
)
from .mosdepth import remove_intermediate_files

from functools import partial
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed


# In[1]: Main Function to Run Normalize Mosdepth
def normalize_mosdepth(config, console):
    """
    Normalizes mosdepth coverage for all samples in the config file.

    Args:
        config (dict): The configuration data loaded from the config file.
        console (Console): The Rich console object for logging.
    """
    try:
        samples_file = config["samples_file"]
        samples = get_samples(samples_file)
        chrom = config.get("chrom", None)
        start = config.get("start_bp", None)
        end = config.get("end_bp", None)
        threads = config.get("threads", 1)
        output_file_prefix = (
            config.get("mosdepth", {}).get("normalize", {}).get("output_file_prefix", None)
        )
        output_file_type = config.get("output_file_type", "tsv")
        output_dir = config.get("output_dir", ".")
        output_file = Path(f"{output_dir}/{output_file_prefix}.{output_file_type}.gz")
        remove_intermediate = config.get("mosdepth", {}).get("remove_intermediate", False)

        threads = config.get("threads", 1)

        mosdepth_dir = config.get("mosdepth", {}).get("work_dir", None)
        min_depth = config["mosdepth"]["normalize"].get("min_depth", 20)
        max_depth = config["mosdepth"]["normalize"].get("max_depth", 100)
        top_frac = config["mosdepth"]["normalize"].get("top_frac", 0.1)
        repeat_mask = config["mosdepth"]["normalize"].get("repeat_mask_file", None)
    except Exception as e:
        log(console, f"[red]Config error: {e}[/red]")
        return

    output_path = Path(output_file).expanduser()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path = setup_output_file(output_path, chrom, start, end)

    # console.rule("[bold blue]Step 1: Collect Individuals")
    individuals = map_mosdepth_files_to_samples(mosdepth_dir, samples)
    if not individuals:
        log(console, f"✗ No mosdepth files found in {mosdepth_dir}", style="danger")
        sys.exit(1)

    # Load repeat mask if provided
    excluded = load_repeat_mask(repeat_mask)

    region_pop_means = compute_population_mean_depths(
        individuals=individuals,
        mosdepth_dir=mosdepth_dir,
        chromosome=chrom,
        start=start,
        end=end,
        excluded=excluded,
        threads=threads,
        console=console,
    )

    valid_regions = {
        region for region, mean_d in region_pop_means.items() if min_depth <= mean_d <= max_depth
    }

    process_func = partial(
        process_one_individual,
        mosdepth_dir=mosdepth_dir,
        chromosome=chrom,
        start=start,
        end=end,
        excluded=excluded,
        valid_regions=valid_regions,
    )

    regions_to_extract = {}
    with progress_bar(
        console, total=len(individuals), description="Extracting per-sample regions..."
    ) as (progress, task):
        with ThreadPoolExecutor(max_workers=max(1, threads)) as executor:
            futures_to_file = {
                executor.submit(process_func, ind_id): ind_id for ind_id in individuals.keys()
            }

            for future in as_completed(futures_to_file):
                ind_id = futures_to_file[future]
                try:
                    ind_id, regions = future.result()
                    regions_to_extract[ind_id] = regions
                except Exception as e:
                    log(console, f"Error processing {ind_id}: {e}", style="danger")
                finally:
                    progress.update(task, advance=1)

    regions_to_extract = filter_empty_samples(regions_to_extract, console)
    if not regions_to_extract:
        log(console, "No valid samples with regions found.", style="danger")
        sys.exit(1)

    individuals_order, mat = build_matrix_from_regions(regions_to_extract)
    individual_raw_means = np.nanmean(mat, axis=1)

    # Normalize the matrix and compute variance ratios
    normalized_mat, variance_ratios, col_means, col_vars = normalize_matrix(mat)

    # keep the top (1 - top_frac) fraction, i.e. everything above the
    # top_frac-th quantile threshold (uses top_frac=0.1 → keeps 90%).
    selected_indices = select_high_variance_regions(variance_ratios, top_frac)

    write_normalized_output(
        normalized_mat,
        individuals_order,
        selected_indices,
        output_path,
        col_means=col_means,
        col_vars=col_vars,
        individual_raw_means=individual_raw_means,
    )
    log(
        console,
        f"Mosdepth normalization complete. Results written to {output_path}",
        style="success",
    )

    if remove_intermediate:
        remove_intermediate_files(mosdepth_dir, console, include_region_bed_gz=True)


def map_mosdepth_files_to_samples(mosdepth_dir: str, samples: list[str]) -> dict[str, Path]:
    """
    Map sample IDs to their mosdepth.global.dist files.

    Args:
        mosdepth_dir: Directory containing mosdepth output files
        samples: List of sample IDs

    Returns:
        Dictionary mapping sample IDs to their global.dist file paths
    """
    mosdepth_dir = Path(mosdepth_dir)
    return {
        f.stem.split(".")[0]: f
        for f in mosdepth_dir.glob("*.regions.bed.gz")
        if f.stem.split(".")[0] in samples
    }


def load_repeat_mask(repeat_bed: str) -> dict[str, set[int]]:
    """
    Load repeat regions into {chrom: set(kb_bins)}.

    Args:
        repeat_bed: Path to repeat mask BED file

    Returns:
        Dictionary mapping chromosomes to sets of excluded kb bins
    """
    excluded = defaultdict(set)
    with open(repeat_bed) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split()
            if len(parts) < 3:
                continue
            chrom, start_str, end_str = parts[0], parts[1], parts[2]

            # FIX #4: normalise chrom to 'chrN' form so matching works
            # regardless of whether the mask uses 'chr6' or '6'.
            chrom = norm_chrom(chrom)

            try:
                start, end = int(start_str), int(end_str)
            except ValueError:
                continue
            for kb in range(start // 1000, end // 1000 + 1):
                excluded[chrom].add(kb)
    return excluded


def norm_chrom(chrom: str) -> str:
    """
    Normalise a chromosome name to 'chrN' form.
    '6' -> 'chr6', 'chr6' -> 'chr6', 'chrX' -> 'chrX'
    """
    return chrom if chrom.startswith("chr") else f"chr{chrom}"


def compute_population_mean_depths(
    individuals: dict[str, Path],
    mosdepth_dir: str,
    chromosome: str,
    start: int,
    end: int,
    excluded: dict[str, set[int]],
    threads: int = 1,
    console=None,
) -> dict[tuple[int, int], float]:
    """
    Compute the population mean depth for every region across all individuals.

    This mirrors the C++ first pass (batches 0–9) that computes mean_depths[]
    and uses it as the depth filter threshold — so the region mask is universal
    (same regions included/excluded for every sample) rather than per-sample.

    Args:
        individuals:  {sample_id: bed_gz_path}
        mosdepth_dir: directory with mosdepth output (passed through to reader)
        chromosome:   target chromosome (or None for all)
        start/end:    target region bounds (or None for whole chromosome)
        excluded:     repeat-mask kb bins {chrom: set(kb)}
        threads:      worker threads for parallel reading
        console:      Rich console for logging

    Returns:
        {(region_start, region_end): population_mean_depth}
    """
    region_sums = defaultdict(float)
    region_counts = defaultdict(int)
    lock_data: dict = {}  # will use a threading.Lock below

    from threading import Lock

    data_lock = Lock()

    def _read_one(ind_id):
        bed_gz = find_bed_gz_for_individual(ind_id, mosdepth_dir)
        if not bed_gz.exists():
            return
        chrom_to_match = norm_chrom(chromosome) if chromosome else None
        local = {}
        try:
            with gzip.open(bed_gz, "rt") as f:
                for line in f:
                    if chrom_to_match and not line.startswith(chrom_to_match):
                        continue
                    fields = line.strip().split("\t")
                    if len(fields) < 4:
                        continue
                    chrom_f = norm_chrom(fields[0])
                    reg_start = int(fields[1])
                    reg_end = int(fields[2])
                    depth = float(fields[3])

                    if start is not None and end is not None:
                        if not (depth > 0 and reg_end >= start and reg_start <= end):
                            continue
                    elif depth <= 0:
                        continue

                    # FIX #4: use normalised chrom when checking repeat mask
                    region_kb = set(range(reg_start // 1000, reg_end // 1000 + 1))
                    if region_kb & excluded.get(chrom_f, set()):
                        continue

                    local[(reg_start, reg_end)] = depth
        except Exception:
            return

        with data_lock:
            for region, d in local.items():
                region_sums[region] += d
                region_counts[region] += 1

    with ThreadPoolExecutor(max_workers=max(1, threads)) as ex:
        list(ex.map(_read_one, individuals.keys()))

    return {
        region: region_sums[region] / region_counts[region]
        for region in region_sums
        if region_counts[region] > 0
    }


def process_one_individual(
    individual_id, mosdepth_dir, chromosome, start, end, valid_regions, excluded
):
    """
    This runs in its own process. Returns (individual_id, [(start,end,depth),...])
    args_tuple contains:

    individual_id, mosdepth_dir, chromosome, start, end, min_depth, max_depth, excluded
    """
    bed_gz = find_bed_gz_for_individual(individual_id, mosdepth_dir)
    if not bed_gz.exists():
        return individual_id, []

    chrom_to_match = norm_chrom(chromosome) if chromosome else None
    results = []

    try:
        with gzip.open(bed_gz, "rt") as f:
            for line in f:
                # Check if line starts with the target chromosome
                if chrom_to_match and not line.startswith(chrom_to_match):
                    continue

                fields = line.strip().split("\t")
                if len(fields) < 4:
                    continue

                chrom = norm_chrom(fields[0])
                region_start = int(fields[1])
                region_end = int(fields[2])
                depth = float(fields[3])

                if start is not None and end is not None:
                    if not (depth > 0 and region_end >= start and region_start <= end):
                        continue
                else:
                    if depth <= 0:
                        continue

                # depth filter
                if (region_start, region_end) not in valid_regions:
                    continue

                # repeat exclusion (kb bins)
                region_kb = set(range(region_start // 1000, region_end // 1000 + 1))
                if region_kb & excluded.get(chrom, set()):
                    continue

                results.append((region_start, region_end, depth))
    except Exception:
        # If file corrupt or other IO error, just return empty
        return individual_id, []

    return individual_id, results


def build_depth_matrix(
    regions_to_extract: dict[str, list[tuple[int, int, float]]],
) -> dict[str, dict[tuple[int, int], float]]:
    """
    Build depth matrix from extracted regions.

    Args:
        regions_to_extract: Dictionary of sample IDs to region lists

    Returns:
        Dictionary mapping {sample: {(start,end): depth}}
    """
    depth_matrix = defaultdict(dict)
    for ind, regions in regions_to_extract.items():
        for start, end, depth in regions:
            depth_matrix[ind][(start, end)] = depth
    return depth_matrix


def build_matrix_from_regions(regions_to_extract, individuals_order=None):
    """
    Build numpy matrix from regions_to_extract.

    Args:
        regions_to_extract: dict of {individual_id: [(start,end,depth),...]}
        individuals_order: optional list of individual IDs to order rows

    Returns:
        individuals_order: list of individual IDs in order
        regions_list: list of (start,end) tuples in order
        mat: numpy array of shape (n_individuals, n_regions) with depths
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

    return individuals_order, mat


def normalize_matrix(mat):
    """
    Normalize the depth matrix to match C++ normalize_mosdepth_inflow logic.

    Steps:
      1) Row-wise (within-individual): divide each row by its mean depth.
      2) Column-wise (across-individual): for each region compute mu and sigma2,
         then transform x -> (x - mu) / sqrt(mu).
      3) Compute variance ratios (100 * sigma2 / mu) for every region.
      4) Median-based rescaling: multiply all values by 1/sqrt(median_ratio/100)
         so the output approximates true z-scores (matches C++ `scale` factor).

    Args:
        mat: numpy array shape (n_individuals, n_regions) with raw depths.

    Returns:
        normalized_mat : numpy array same shape, rescaled z-scores.
        variance_ratios: dict {region_index: 100*sigma2/mu} for non-NaN regions.
    """
    mat = mat.copy()

    row_means = np.nanmean(mat, axis=1)
    row_means_safe = np.where(row_means == 0, np.nan, row_means)
    mat = (mat.T / row_means_safe).T

    n_inds = mat.shape[0]
    col_means = np.nanmean(mat, axis=0)  # mu
    col_vars = np.nansum((mat - col_means) ** 2, axis=0) / (n_inds - 1)  # s2, ddof=1

    # variance ratio: ratioMult * s2 / mu
    ratio_mult = 100.0
    with np.errstate(invalid="ignore", divide="ignore"):
        var_ratio = np.where(col_means > 0, ratio_mult * col_vars / col_means, np.nan)

    # transform: (x - mu) / sqrt(mu)
    mu_pos = col_means > 0
    sqrt_mu = np.where(
        mu_pos, np.sqrt(col_means, where=mu_pos, out=np.full_like(col_means, np.nan)), np.nan
    )
    mat[:, mu_pos] = (mat[:, mu_pos] - col_means[mu_pos]) / sqrt_mu[mu_pos]

    valid_ratios = var_ratio[~np.isnan(var_ratio)]
    if valid_ratios.size > 0:
        sigma2ratio_median = float(np.median(valid_ratios))
        if sigma2ratio_median > 0:
            scale = 1.0 / np.sqrt(sigma2ratio_median / ratio_mult)
        else:
            scale = 1.0
    else:
        scale = 1.0

    mat *= scale  # apply rescaling to the whole matrix

    variance_ratios = {
        i: float(var_ratio[i]) for i in range(len(col_means)) if not np.isnan(var_ratio[i])
    }

    return mat, variance_ratios, col_means, col_vars  # expose means/vars for header


def select_high_variance_regions(
    variance_ratios: dict[int, float], top_frac: float = 0.9
) -> list[int]:
    """
    Select top fraction of regions by variance ratio.

    Args:
        variance_ratios: {region_index: variance_ratio}
        top_frac: fraction of top regions to retain

    Returns:
        List of selected region indices
    """
    if not variance_ratios:
        return []

    sorted_ratios = sorted(variance_ratios.values())
    threshold_idx = int(top_frac * len(sorted_ratios))
    threshold = sorted_ratios[threshold_idx]

    return [idx for idx, ratio in variance_ratios.items() if ratio > threshold]


def write_normalized_output(
    mat: np.ndarray,
    individuals_order: list,
    selected_indices: list,
    output_file: Path,
    col_means: np.ndarray,
    col_vars: np.ndarray,
    individual_raw_means: np.ndarray,
    ratio_mult: float = 100.0,
):
    """
    Write normalised matrix in the C++-compatible format.

    File layout
    -----------
    Line 0  : N <tab> Rwant <tab> mu_1 <tab> mu_2 ...
    Line 1  : N <tab> Rwant <tab> varRatio_1 <tab> varRatio_2 ...
    Line 2+ : ID <tab> indiv_scale <tab> z_1 <tab> z_2 ...

    Args:
        mat                  : normalised + rescaled matrix (n_individuals × n_regions).
        individuals_order    : list of sample IDs.
        regions_list         : list of (start, end) tuples for all regions.
        selected_indices     : indices of high-variance regions to keep.
        output_file          : Path for the output .tsv.gz file.
        col_means            : per-region means (length = n_regions).
        col_vars             : per-region variances ddof=1 (length = n_regions).
        individual_raw_means : per-sample mean depth before any normalisation
                               (used to write the `scale` column; matches C++
                               `0.01f*scales[n]` where scales[n] was in units
                               of 0.01x, so the written value is in 1x units).
        ratio_mult           : multiplier used for variance ratios (default 100).
    """
    N = len(individuals_order)
    Rwant = len(selected_indices)

    sel_means = col_means[selected_indices]
    sel_vars = col_vars[selected_indices]
    with np.errstate(invalid="ignore", divide="ignore"):
        sel_ratios = np.where(sel_means > 0, ratio_mult * sel_vars / sel_means, np.nan)

    with gzip.open(output_file, "wt") as out:
        means_str = "\t".join("NA" if np.isnan(v) else f"{v:.3f}" for v in sel_means)
        out.write(f"{N}\t{Rwant}\t{means_str}\n")

        ratios_str = "\t".join("NA" if np.isnan(v) else f"{v:.3f}" for v in sel_ratios)
        out.write(f"{N}\t{Rwant}\t{ratios_str}\n")

        for i, ind_id in enumerate(individuals_order):
            # FIX #3: write raw mean directly (already 1x, no *0.01 needed)
            ind_scale = individual_raw_means[i]
            vals = ["NA" if np.isnan(mat[i, j]) else f"{mat[i, j]:.2f}" for j in selected_indices]
            out.write(f"{ind_id}\t{ind_scale:.2f}\t" + "\t".join(vals) + "\n")


def find_bed_gz_for_individual(individual_id: str, mosdepth_dir: str) -> Path:
    """
    Find the .regions.bed.gz file for a given individual in the mosdepth directory.

    Args:
        individual_id: Sample ID to find
        mosdepth_dir: Directory containing mosdepth output files

    Returns:
        Path to the individual's .regions.bed.gz file, or Path to non-existent file if not found
    """
    mosdepth_dir = Path(mosdepth_dir)
    matches = list(mosdepth_dir.glob(f"*{individual_id}*regions.bed.gz"))
    if matches:
        return matches[0]
    else:
        return mosdepth_dir / f"{individual_id}.regions.bed.gz"


def filter_empty_samples(regions_to_extract: dict[str, list[tuple[int, int, float]]], console=None):
    """
    Remove samples that have zero regions.

    Args:
        regions_to_extract: dict of {sample_id: [(start, end, depth), ...]}
        console: optional Rich console for logging

    Returns:
        Filtered dictionary with only samples that have ≥1 region
    """
    n_before = len(regions_to_extract)

    filtered = {
        sample: regions for sample, regions in regions_to_extract.items() if len(regions) > 0
    }

    n_after = len(filtered)
    n_removed = n_before - n_after

    if n_removed > 0:
        msg = f"Removed {n_removed} samples with 0 regions"
        log(console, msg, style="warning") if console else print(msg)

    return filtered
