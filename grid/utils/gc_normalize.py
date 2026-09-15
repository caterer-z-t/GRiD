"""
GC-content LOWESS normalization for KIV-2 copy number estimation.

Implements the key steps from the DRAGEN LPA caller (Behera et al. 2024):
  1. For each sample, read per-region depths from mosdepth output.
  2. Pair depths with GC content from a precomputed normalization-regions file.
  3. Fit a per-sample LOWESS curve: depth_per_bp ~ GC_content.
  4. Predict expected depth at the KIV-2 region's GC content.
  5. GC-corrected depth ratio = observed_KIV2_depth / predicted_depth.
  6. Scale by reference_cn (default 6.2, following DRAGEN's empirical calibration)
     to produce total KIV-2 copy number.
  7. Compute MAD-based QC score across all regions for flagging.

Output TSV columns:
  sample_id  gc_corrected_cn  diploid_cn  mad_score  qc_flag
where
  gc_corrected_cn  = GC-LOWESS-corrected total copy number
  diploid_cn       = same (alias retained for pipeline compatibility)
  mad_score        = median absolute deviation of normalised depths (DRAGEN QC)
  qc_flag          = 1 if mad_score > mad_threshold else 0
"""
from __future__ import annotations

import gzip
import math
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from functools import partial
from pathlib import Path
from threading import Lock

import numpy as np

from .utils import log, get_samples, progress_bar
from .compute_gc_content import load_gc_content_file, compute_gc_for_regions
from .normalize_mosdepth import norm_chrom


# ---------------------------------------------------------------------------
# Pipeline entry point
# ---------------------------------------------------------------------------

def gc_normalize(config: dict, console=None) -> None:
    """
    Pipeline step: GC-LOWESS normalization for all samples.

    Reads:
      - mosdepth .regions.bed.gz files  (mosdepth.work_dir)
      - normalization regions GC file   (gc_normalize.norm_regions_gc_file)

    Writes:
      - TSV with per-sample GC-corrected CN  (gc_normalize.output_file_prefix)
    """
    try:
        gc_cfg = config.get("gc_normalize", {})
        output_dir = config.get("output_dir", ".")
        output_prefix = gc_cfg.get("output_file_prefix", "gc_corrected_cn")
        output_file_type = config.get("output_file_type", "tsv")
        output_file = Path(output_dir) / f"{output_prefix}.{output_file_type}"

        mosdepth_dir = config.get("mosdepth", {}).get("work_dir")
        if not mosdepth_dir:
            log(console, "gc_normalize requires mosdepth.work_dir in config", style="danger")
            sys.exit(1)

        norm_regions_gc_file = gc_cfg.get("norm_regions_gc_file")
        if not norm_regions_gc_file:
            log(
                console,
                "gc_normalize requires gc_normalize.norm_regions_gc_file in config.\n"
                "Run the select_norm_regions step first to generate this file.",
                style="danger",
            )
            sys.exit(1)

        target_chrom = config.get("chrom")
        target_start = config.get("start_bp")
        target_end = config.get("end_bp")
        reference_genome = config.get("reference_genome")

        reference_cn = gc_cfg.get("reference_cn", 6.2)
        lowess_frac = gc_cfg.get("lowess_frac", 0.2)
        mad_threshold = gc_cfg.get("mad_threshold", 0.11)
        threads = config.get("threads", 1)
        samples_file = config["samples_file"]
    except Exception as e:
        log(console, f"Config error: {e}", style="danger")
        sys.exit(1)

    output_file.parent.mkdir(parents=True, exist_ok=True)

    # --- Load normalization regions and their GC content ---
    log(console, f"Loading normalization regions from {norm_regions_gc_file}…", style="info")
    gc_map = load_gc_content_file(norm_regions_gc_file)
    if not gc_map:
        log(console, f"No regions loaded from {norm_regions_gc_file}", style="danger")
        sys.exit(1)
    log(console, f"Loaded GC content for {len(gc_map)} normalization regions.", style="info")

    # --- Compute KIV-2 region GC content ---
    target_chrom_norm = norm_chrom(target_chrom) if target_chrom else None
    kiv2_gc = _compute_kiv2_gc(
        reference_genome, target_chrom, target_start, target_end, console
    )
    if math.isnan(kiv2_gc):
        log(
            console,
            "Could not compute GC content for target region; "
            "check reference_genome and chrom/start_bp/end_bp in config.",
            style="danger",
        )
        sys.exit(1)
    log(
        console,
        f"Target region GC content: {kiv2_gc:.3f} "
        f"({target_chrom}:{target_start}-{target_end})",
        style="info",
    )

    # --- Discover per-sample mosdepth files ---
    samples = get_samples(samples_file)
    mosdepth_dir_path = Path(mosdepth_dir)
    sample_files = _map_samples_to_mosdepth_files(mosdepth_dir_path, samples)
    log(console, f"Found mosdepth files for {len(sample_files)}/{len(samples)} samples.", style="info")

    # --- Process each sample ---
    process_fn = partial(
        _process_one_sample,
        gc_map=gc_map,
        target_chrom=target_chrom_norm,
        target_start=target_start,
        target_end=target_end,
        kiv2_gc=kiv2_gc,
        reference_cn=reference_cn,
        lowess_frac=lowess_frac,
        mad_threshold=mad_threshold,
    )

    results: list[tuple[str, float, float, float, int]] = []
    write_lock = Lock()

    with progress_bar(console, total=len(sample_files), description="GC normalising samples…") as (
        progress,
        task,
    ):
        with ThreadPoolExecutor(max_workers=max(1, threads)) as executor:
            future_to_sample = {
                executor.submit(process_fn, bed_gz): sample_id
                for sample_id, bed_gz in sample_files.items()
            }
            for future in as_completed(future_to_sample):
                sample_id = future_to_sample[future]
                try:
                    result = future.result()
                    if result is not None:
                        results.append((sample_id, *result))
                except Exception as exc:
                    log(console, f"Error processing {sample_id}: {exc}", style="warning")
                finally:
                    progress.update(task, advance=1)

    # --- Write output ---
    results.sort(key=lambda r: r[0])
    output_file.parent.mkdir(parents=True, exist_ok=True)
    with open(output_file, "w") as out:
        out.write("sample_id\tgc_corrected_cn\tdiploid_cn\tmad_score\tqc_flag\n")
        for sample_id, cn, mad_score, qc_flag in results:
            out.write(f"{sample_id}\t{cn:.4f}\t{cn:.4f}\t{mad_score:.4f}\t{qc_flag}\n")

    n_flagged = sum(1 for _, _cn, mad, _flag in results if mad > mad_threshold)
    log(
        console,
        f"GC normalisation complete. {len(results)} samples processed, "
        f"{n_flagged} flagged (MAD > {mad_threshold}).\n"
        f"Results written to {output_file}",
        style="success",
    )


# ---------------------------------------------------------------------------
# Per-sample processing
# ---------------------------------------------------------------------------

def _process_one_sample(
    bed_gz: Path,
    gc_map: dict[tuple[str, int, int], float],
    target_chrom: str | None,
    target_start: int | None,
    target_end: int | None,
    kiv2_gc: float,
    reference_cn: float,
    lowess_frac: float,
    mad_threshold: float,
) -> tuple[float, float, int] | None:
    """
    For one sample: fit LOWESS, correct KIV-2 depth, return (cn, mad_score, qc_flag).
    Returns None on failure.
    """
    norm_depths, kiv2_depth = _read_depths_from_mosdepth(
        bed_gz, gc_map, target_chrom, target_start, target_end
    )
    if len(norm_depths) < 10:
        return None
    if kiv2_depth is None or kiv2_depth <= 0:
        return None

    gcs = np.array([gc for gc, _ in norm_depths])
    depths = np.array([d for _, d in norm_depths])

    # Fit the LOWESS grid once; reuse for KIV-2 prediction and MAD QC
    try:
        grid_x, grid_y = _lowess_fit(gcs, depths, lowess_frac)
    except Exception:
        return None

    predicted_depth = float(np.interp(kiv2_gc, grid_x, grid_y))
    if not np.isfinite(predicted_depth) or predicted_depth <= 0:
        return None

    # GC-corrected depth ratio → total CN
    cn = (kiv2_depth / predicted_depth) * reference_cn

    # QC: MAD of (observed / expected) across all normalization regions.
    # Expected depth at each region's GC content comes from the same LOWESS fit.
    predicted_at_norm = np.interp(gcs, grid_x, grid_y)
    with np.errstate(divide="ignore", invalid="ignore"):
        norm_factors = np.where(predicted_at_norm > 0, depths / predicted_at_norm, np.nan)
    norm_factors = norm_factors[np.isfinite(norm_factors)]
    mad_score = _mad(norm_factors) if norm_factors.size > 0 else math.nan
    qc_flag = 1 if (np.isfinite(mad_score) and mad_score > mad_threshold) else 0

    return cn, mad_score, qc_flag


def _read_depths_from_mosdepth(
    bed_gz: Path,
    gc_map: dict[tuple[str, int, int], float],
    target_chrom: str | None,
    target_start: int | None,
    target_end: int | None,
) -> tuple[list[tuple[float, float]], float | None]:
    """
    Read mosdepth regions.bed.gz for one sample.

    Returns:
        norm_depths : list of (gc_fraction, per-bp_depth) for normalization regions
        kiv2_depth  : weighted mean per-bp depth over the KIV-2 target region, or None
    """
    norm_depths: list[tuple[float, float]] = []
    kiv2_cov = 0.0
    kiv2_bp = 0

    try:
        with gzip.open(bed_gz, "rt") as f:
            for line in f:
                fields = line.strip().split("\t")
                if len(fields) < 4:
                    continue
                chrom = norm_chrom(fields[0])
                start = int(fields[1])
                end = int(fields[2])
                depth = float(fields[3])
                region_len = end - start
                if region_len <= 0:
                    continue

                per_bp = depth  # mosdepth reports mean per-base depth already

                # Check if this bin overlaps the KIV-2 target region
                if (
                    target_chrom
                    and chrom == target_chrom
                    and target_start is not None
                    and target_end is not None
                    and not (end <= target_start or start >= target_end)
                ):
                    overlap_start = max(start, target_start)
                    overlap_end = min(end, target_end)
                    overlap = overlap_end - overlap_start
                    if overlap > 0:
                        kiv2_cov += per_bp * overlap
                        kiv2_bp += overlap
                    continue

                # Check if this region is in our normalization set
                key = (chrom, start, end)
                if key in gc_map:
                    gc = gc_map[key]
                    if not math.isnan(gc) and per_bp > 0:
                        norm_depths.append((gc, per_bp))
    except Exception:
        return [], None

    kiv2_depth = (kiv2_cov / kiv2_bp) if kiv2_bp > 0 else None
    return norm_depths, kiv2_depth


def _lowess_fit(
    gcs: np.ndarray,
    depths: np.ndarray,
    frac: float,
    n_grid: int = 100,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Fit a LOWESS curve using a pure-NumPy grid approach.

    Fits weighted local linear regression at `n_grid` evenly spaced GC values
    spanning the observed range, then returns the (grid_x, grid_y) arrays for
    subsequent np.interp calls.  This avoids the statsmodels / binary-compat
    issues and runs in O(n_grid × k) time instead of O(n × k).
    """
    n = len(gcs)
    k = max(2, int(frac * n))
    gc_min, gc_max = float(np.min(gcs)), float(np.max(gcs))
    if gc_max == gc_min:
        return np.array([gc_min, gc_max]), np.array([float(np.mean(depths))] * 2)

    grid_x = np.linspace(gc_min, gc_max, n_grid)
    grid_y = np.empty(n_grid)

    for i, x0 in enumerate(grid_x):
        dists = np.abs(gcs - x0)
        idx = np.argpartition(dists, k - 1)[:k]  # k nearest (unsorted)
        h = float(np.max(dists[idx]))
        if h < 1e-10:
            grid_y[i] = float(np.mean(depths[idx]))
            continue
        u = dists[idx] / h
        w = np.clip((1.0 - u ** 3) ** 3, 0.0, None)
        x_nbr = gcs[idx]
        y_nbr = depths[idx]
        sw = np.sum(w)
        swx = np.dot(w, x_nbr)
        swx2 = np.dot(w, x_nbr ** 2)
        swy = np.dot(w, y_nbr)
        swxy = np.dot(w, x_nbr * y_nbr)
        det = sw * swx2 - swx * swx
        if abs(det) < 1e-12:
            grid_y[i] = float(np.dot(w, y_nbr) / sw) if sw > 0 else float(np.mean(y_nbr))
        else:
            beta0 = (swx2 * swy - swx * swxy) / det
            beta1 = (sw * swxy - swx * swy) / det
            grid_y[i] = float(beta0 + beta1 * x0)

    return grid_x, grid_y


def _lowess_predict(
    gcs: np.ndarray,
    depths: np.ndarray,
    target_gc: float,
    frac: float,
    n_grid: int = 100,
) -> float | None:
    """
    Fit LOWESS on (gcs, depths) and return the predicted depth at target_gc.
    Returns None if fitting fails or prediction is non-finite.
    """
    try:
        grid_x, grid_y = _lowess_fit(gcs, depths, frac, n_grid)
        predicted = float(np.interp(target_gc, grid_x, grid_y))
        return predicted if np.isfinite(predicted) else None
    except Exception:
        return None


def _mad(arr: np.ndarray) -> float:
    """Median absolute deviation."""
    if arr.size == 0:
        return math.nan
    return float(np.median(np.abs(arr - np.median(arr))))


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _compute_kiv2_gc(
    reference_genome: str | None,
    chrom: str | None,
    start: int | None,
    end: int | None,
    console=None,
) -> float:
    """Compute GC content for the KIV-2 target region from the reference FASTA."""
    if not reference_genome or not chrom or start is None or end is None:
        return math.nan
    try:
        import pysam
        with pysam.FastaFile(str(reference_genome)) as fasta:
            from .compute_gc_content import compute_gc_for_region
            return compute_gc_for_region(fasta, chrom, start, end)
    except Exception as e:
        log(console, f"Failed to compute KIV-2 GC content: {e}", style="warning")
        return math.nan


def _map_samples_to_mosdepth_files(
    mosdepth_dir: Path, samples: list[str]
) -> dict[str, Path]:
    """Map sample IDs to their mosdepth .regions.bed.gz files."""
    sample_set = set(samples)
    result: dict[str, Path] = {}
    for bed_gz in mosdepth_dir.glob("*.regions.bed.gz"):
        name_part = bed_gz.name.split(".")[0]
        parts = name_part.split("_")
        for i in range(len(parts), 0, -1):
            candidate = "_".join(parts[:i])
            if candidate in sample_set and candidate not in result:
                result[candidate] = bed_gz
                break
    return result
