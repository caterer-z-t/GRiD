"""
Select stable, GC-diverse normalization regions for GC LOWESS correction.

Strategy (mirrors DRAGEN's approach):
  1. Read per-region depths across all samples from existing mosdepth output.
  2. Compute per-region population mean and CV.
  3. Filter: pop_mean in [min_depth, max_depth] AND CV < max_cv.
  4. Exclude the target VNTR region and any user-supplied blacklist regions.
  5. Compute GC content for all passing regions from the reference FASTA.
  6. Bin by GC content and select the top `n_per_gc_bin` lowest-CV regions per bin.
  7. Write selected regions (chrom, start, end, gc_fraction) to output file.

This step is intended to run once per cohort and its output is then referenced
by the gc_normalize step for every subsequent sample.
"""
from __future__ import annotations

import gzip
import sys
import math
from collections import defaultdict
from pathlib import Path
from threading import Lock
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np

from .utils import log, get_samples, progress_bar
from .compute_gc_content import compute_gc_for_regions, write_gc_content_file
from .normalize_mosdepth import norm_chrom


def select_normalization_regions(config: dict, console=None) -> None:
    """
    Pipeline entry point: select stable, GC-diverse normalization regions.

    Writes a TSV file with columns: chrom, start, end, gc_fraction.
    This file is consumed by the gc_normalize step.
    """
    try:
        snr_cfg = config.get("select_norm_regions", {})
        output_dir = config.get("output_dir", ".")
        output_file = Path(output_dir) / snr_cfg.get("output_file", "norm_regions.gc.tsv")

        mosdepth_dir = config.get("mosdepth", {}).get("work_dir")
        if not mosdepth_dir:
            log(console, "select_norm_regions requires mosdepth.work_dir", style="danger")
            sys.exit(1)

        reference_genome = config.get("reference_genome")
        if not reference_genome:
            log(console, "select_norm_regions requires reference_genome", style="danger")
            sys.exit(1)

        samples_file = config["samples_file"]
        target_chrom = config.get("chrom")
        target_start = config.get("start_bp")
        target_end = config.get("end_bp")

        min_depth = config.get("mosdepth", {}).get("normalize", {}).get("min_depth", 20)
        max_depth = config.get("mosdepth", {}).get("normalize", {}).get("max_depth", 100)
        max_cv = snr_cfg.get("max_cv", 0.15)
        gc_bin_size = snr_cfg.get("gc_bin_size", 0.05)
        n_per_gc_bin = snr_cfg.get("n_per_gc_bin", 300)
        threads = config.get("threads", 1)
    except Exception as e:
        log(console, f"Config error: {e}", style="danger")
        sys.exit(1)

    output_file.parent.mkdir(parents=True, exist_ok=True)

    from .utils import get_samples
    samples = get_samples(samples_file)

    # --- Step 1: collect per-region depth sums / sum-of-squares across samples ---
    log(console, f"Reading mosdepth files from {mosdepth_dir}…", style="info")
    region_sums, region_sq_sums, region_counts = _collect_region_stats(
        mosdepth_dir=mosdepth_dir,
        samples=samples,
        target_chrom=norm_chrom(target_chrom) if target_chrom else None,
        target_start=target_start,
        target_end=target_end,
        min_depth=min_depth,
        max_depth=max_depth,
        threads=threads,
        console=console,
    )

    if not region_sums:
        log(console, "No usable regions found in mosdepth output.", style="danger")
        sys.exit(1)

    # --- Step 2: compute mean and CV, filter ---
    log(console, "Computing per-region statistics and filtering…", style="info")
    passing_regions = _filter_regions(
        region_sums, region_sq_sums, region_counts,
        min_depth=min_depth, max_depth=max_depth, max_cv=max_cv,
    )

    log(
        console,
        f"{len(passing_regions)} regions passed depth+CV filters "
        f"(from {len(region_sums)} total).",
        style="info",
    )

    if not passing_regions:
        log(console, "No regions survived filters; cannot select normalization regions.", style="danger")
        sys.exit(1)

    # --- Step 3: compute GC content for passing regions ---
    log(console, f"Computing GC content from {reference_genome}…", style="info")
    gc_map = compute_gc_for_regions(reference_genome, [r for r, _ in passing_regions])

    # --- Step 4: bin by GC and select top n_per_gc_bin lowest-CV per bin ---
    selected = _select_gc_diverse(passing_regions, gc_map, gc_bin_size, n_per_gc_bin)

    log(
        console,
        f"Selected {len(selected)} normalization regions covering "
        f"{len({int(gc_map.get(r, 0) / gc_bin_size) for r in selected})} GC bins.",
        style="success",
    )

    # --- Step 5: write output ---
    selected_gc = {r: gc_map[r] for r in selected if r in gc_map}
    write_gc_content_file(selected_gc, output_file)
    log(console, f"Normalization regions written to {output_file}", style="success")


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _collect_region_stats(
    mosdepth_dir: str,
    samples: list[str],
    target_chrom: str | None,
    target_start: int | None,
    target_end: int | None,
    min_depth: float,
    max_depth: float,
    threads: int,
    console=None,
) -> tuple[dict, dict, dict]:
    """
    First pass over mosdepth files: accumulate sum, sum-of-squares, and count
    per region so we can compute mean and std in one pass.

    Regions that overlap the target VNTR are excluded.
    """
    region_sums: dict[tuple[str, int, int], float] = defaultdict(float)
    region_sq_sums: dict[tuple[str, int, int], float] = defaultdict(float)
    region_counts: dict[tuple[str, int, int], int] = defaultdict(int)
    lock = Lock()

    mosdepth_dir_path = Path(mosdepth_dir)

    def _read_one(sample_id: str) -> None:
        matches = list(mosdepth_dir_path.glob(f"*{sample_id}*regions.bed.gz"))
        if not matches:
            return
        bed_gz = matches[0]
        local_sums: dict = {}
        local_sq: dict = {}
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

                    if depth <= 0:
                        continue

                    # Exclude the VNTR target region
                    if (
                        target_chrom
                        and chrom == target_chrom
                        and target_start is not None
                        and target_end is not None
                        and not (end <= target_start or start >= target_end)
                    ):
                        continue

                    key = (chrom, start, end)
                    local_sums[key] = local_sums.get(key, 0.0) + depth
                    local_sq[key] = local_sq.get(key, 0.0) + depth * depth
        except Exception:
            return

        with lock:
            for key, s in local_sums.items():
                region_sums[key] += s
                region_sq_sums[key] += local_sq[key]
                region_counts[key] += 1

    with ThreadPoolExecutor(max_workers=max(1, threads)) as ex:
        list(ex.map(_read_one, samples))

    return dict(region_sums), dict(region_sq_sums), dict(region_counts)


def _filter_regions(
    sums: dict,
    sq_sums: dict,
    counts: dict,
    min_depth: float,
    max_depth: float,
    max_cv: float,
) -> list[tuple[tuple[str, int, int], float]]:
    """
    Return list of (region, cv) for regions passing depth and CV filters.
    """
    passing = []
    for key in sums:
        n = counts.get(key, 0)
        if n < 2:
            continue
        mean = sums[key] / n
        if not (min_depth <= mean <= max_depth):
            continue
        # Unbiased sample variance
        var = (sq_sums[key] - n * mean * mean) / (n - 1)
        if var < 0:
            var = 0.0
        std = math.sqrt(var)
        cv = std / mean if mean > 0 else math.inf
        if cv <= max_cv:
            passing.append((key, cv))
    return passing


def _select_gc_diverse(
    passing_regions: list[tuple[tuple[str, int, int], float]],
    gc_map: dict[tuple[str, int, int], float],
    gc_bin_size: float,
    n_per_gc_bin: int,
) -> list[tuple[str, int, int]]:
    """
    Bin regions by GC content and select up to n_per_gc_bin lowest-CV regions
    per bin.  Regions without a valid GC value are skipped.
    """
    bins: dict[int, list[tuple[float, tuple[str, int, int]]]] = defaultdict(list)
    for region, cv in passing_regions:
        gc = gc_map.get(region, math.nan)
        if math.isnan(gc):
            continue
        bin_idx = int(gc / gc_bin_size)
        bins[bin_idx].append((cv, region))

    selected = []
    for bin_idx, entries in bins.items():
        entries.sort(key=lambda x: x[0])  # ascending CV → lowest first
        selected.extend(r for _, r in entries[:n_per_gc_bin])

    return selected
