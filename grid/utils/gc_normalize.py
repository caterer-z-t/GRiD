# In[0]: Imports
from __future__ import annotations

from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from threading import Lock
from functools import partial
import numpy as np
import pandas as pd
import pysam
from statsmodels.nonparametric.smoothers_lowess import lowess

from .utils import log, get_samples, find_file, progress_bar


# In[1]: Main function
def gc_normalize(config, console=None):
    """
    LOWESS GC normalization of LPA read counts across all samples.

    For each sample:
      1. Count reads in background windows spanning the region of interest
      2. Compute RPK (reads per kilobase) for each window and the target region
      3. Fit a LOWESS curve (RPK ~ GC content) to background windows
      4. Compute corrected_rpk_ratio = observed_RPK / predicted_RPK for target
      5. Write per-sample results; aggregate into a combined TSV

    The corrected_rpk_ratio is the value used downstream in compute_diploid_genotypes
    as a GC-corrected substitute for raw read counts.

    Nested under mosdepth in config:
      mosdepth:
        gc_normalize:
          run: True
          output_file_prefix: "gc_normalized"

          # Target region (defaults to top-level chrom/start_bp/end_bp)
          extract_chrom: "chr6"
          extract_start: 160605062
          extract_end:   160640362
          extract_length: 35300       # bp used for RPK denominator

          # KIV-2 reference scaling (DRAGEN method)
          kiv2_ref_length: 33281      # bp spanned by KIV-2 in reference
          kiv2_scale_factor: 6.2      # 6 copies + 0.2 empirical adjustment

          # Background windows
          bg_span_start: 135600000
          bg_span_end:   185600000
          bg_excl_start: 160100000
          bg_excl_end:   161100000
          bg_window_size: 1000
          bg_num_windows: 200

          # LOWESS parameters
          lowess_frac: 0.3
          min_mapq: 0                 # keep multi-mappers for VNTR signal

          # MAD QC threshold (from DRAGEN paper)
          mad_threshold: 0.11
    """
    try:
        directory_loc = config["directory_loc"]
        samples_file = config["samples_file"]
        samples = get_samples(samples_file)
        threads = config.get("threads", 1)
        ref = config.get("reference_genome")

        gc_cfg = config.get("mosdepth", {}).get("gc_normalize", {})
        output_file_prefix = gc_cfg.get("output_file_prefix", "gc_normalized")
        output_file_type = config.get("output_file_type", "tsv")
        output_dir = Path(config.get("output_dir", "."))
        output_file = output_dir / f"{output_file_prefix}.{output_file_type}"

        # Target region — default to top-level chrom/start_bp/end_bp
        extract_chrom = gc_cfg.get("extract_chrom", config.get("chrom", "chr6"))
        extract_start = int(gc_cfg.get("extract_start", config.get("start_bp", 160605062)))
        extract_end = int(gc_cfg.get("extract_end", config.get("end_bp", 160640362)))
        extract_length = int(gc_cfg.get("extract_length", extract_end - extract_start))

        # KIV-2 scaling parameters (DRAGEN method)
        kiv2_ref_length = int(gc_cfg.get("kiv2_ref_length", 33281))
        kiv2_scale_factor = float(gc_cfg.get("kiv2_scale_factor", 6.2))

        # Background window parameters
        bg_span_start = int(gc_cfg.get("bg_span_start", 135_600_000))
        bg_span_end = int(gc_cfg.get("bg_span_end", 185_600_000))
        bg_excl_start = int(gc_cfg.get("bg_excl_start", 160_100_000))
        bg_excl_end = int(gc_cfg.get("bg_excl_end", 161_100_000))
        bg_window_size = int(gc_cfg.get("bg_window_size", 1000))
        bg_num_windows = int(gc_cfg.get("bg_num_windows", 200))

        # LOWESS parameters
        lowess_frac = float(gc_cfg.get("lowess_frac", 0.3))
        min_mapq = int(gc_cfg.get("min_mapq", 0))
        mad_threshold = float(gc_cfg.get("mad_threshold", 0.11))

    except Exception as e:
        log(console, f"Config error: {e}", style="danger")
        return

    output_dir.mkdir(parents=True, exist_ok=True)
    output_file.unlink(missing_ok=True)

    files = {
        sample: result
        for sample in samples
        if (result := find_file(directory_loc, sample, config.get("file_type"))) is not None
    }

    if not files:
        log(console, "No CRAM files found — check directory_loc and samples_file.", style="danger")
        return

    # Pre-compute background regions and GC content once (shared across all samples)
    log(console, "Computing background region coordinates and GC content...", style="info")
    bg_regions = _get_background_regions(
        bg_span_start,
        bg_span_end,
        bg_excl_start,
        bg_excl_end,
        bg_window_size,
        bg_num_windows,
    )
    bg_gc = _compute_gc_for_regions(bg_regions + [(extract_chrom, extract_start, extract_end)], ref)
    extract_gc = bg_gc.pop((extract_chrom, extract_start, extract_end), None)

    if extract_gc is None:
        log(console, "Failed to compute GC content for target region.", style="danger")
        return

    log(
        console,
        f"Target region GC: {extract_gc:.2f}% | Background windows: {len(bg_regions)}",
        style="info",
    )

    write_lock = Lock()
    header_written = [False]

    process_func = partial(
        _process_single_cram,
        ref_fasta=ref,
        extract_chrom=extract_chrom,
        extract_start=extract_start,
        extract_end=extract_end,
        extract_length=extract_length,
        extract_gc=extract_gc,
        bg_regions=bg_regions,
        bg_gc=bg_gc,
        kiv2_ref_length=kiv2_ref_length,
        kiv2_scale_factor=kiv2_scale_factor,
        lowess_frac=lowess_frac,
        min_mapq=min_mapq,
        mad_threshold=mad_threshold,
        console=console,
    )

    failed = []
    with progress_bar(console, total=len(files), description="GC normalizing") as (progress, task):
        with ThreadPoolExecutor(max_workers=threads) as executor:
            future_to_sample = {
                executor.submit(process_func, cram_path): sample
                for sample, cram_path in files.items()
            }
            for future in as_completed(future_to_sample):
                sample = future_to_sample[future]
                result = future.result()
                if result is not None:
                    _write_result(output_file, result, write_lock, header_written)
                else:
                    failed.append(sample)
                    log(console, f"✗ {sample} failed", style="danger")
                progress.advance(task)

    if failed:
        log(console, f"{len(failed)} samples failed: {failed}", style="warning")

    log(console, f"GC normalization complete. Output: {output_file}", style="success")


# In[2]: Process a single CRAM
def _process_single_cram(
    cram_path: str,
    ref_fasta: str,
    extract_chrom: str,
    extract_start: int,
    extract_end: int,
    extract_length: int,
    extract_gc: float,
    bg_regions: list[tuple[str, int, int]],
    bg_gc: dict[tuple[str, int, int], float],
    kiv2_ref_length: int,
    kiv2_scale_factor: float,
    lowess_frac: float,
    min_mapq: int,
    mad_threshold: float,
    console=None,
) -> dict | None:
    """
    Run the full LOWESS GC normalization pipeline for one sample.
    Returns a dict of results or None on failure.
    """
    sample_id = Path(cram_path).stem
    try:
        # Count reads in background windows
        bg_records = []
        for region in bg_regions:
            chrom, start, end = region
            rc = _count_reads(cram_path, ref_fasta, chrom, start, end, min_mapq)
            length = end - start
            bg_records.append(
                {
                    "gc": bg_gc.get(region, float("nan")),
                    "rpk": _rpk(rc, length),
                }
            )
        bg_df = pd.DataFrame(bg_records).dropna()

        # Count reads in target region
        lpa_rc = _count_reads(
            cram_path, ref_fasta, extract_chrom, extract_start, extract_end, min_mapq
        )
        lpa_rpk = _rpk(lpa_rc, extract_length)

        # LOWESS GC correction
        lowess_result = _lowess_gc_correct(extract_gc, lpa_rpk, bg_df, lowess_frac)

        # KIV-2 copy number scaling (DRAGEN method)
        cn_result = _compute_copy_number(
            lowess_result["corrected_rpk"],
            bg_df["rpk"].mean(),
            extract_length,
            kiv2_ref_length,
            kiv2_scale_factor,
        )

        # MAD QC — flag samples with high coverage variation
        all_rpks = bg_df["rpk"].values
        mad = float(np.median(np.abs(all_rpks - np.median(all_rpks))))
        mad_flag = mad > mad_threshold

        return {
            "subject_id": sample_id,
            "lpa_read_count": lpa_rc,
            "lpa_rpk": round(lpa_rpk, 4),
            "lpa_gc": round(extract_gc, 4),
            "predicted_rpk": round(lowess_result["predicted_rpk"], 4),
            "corrected_rpk_ratio": round(lowess_result["corrected_rpk"], 4),
            "bg_residual_mean": round(lowess_result["bg_residual_mean"], 4),
            "bg_residual_std": round(lowess_result["bg_residual_std"], 4),
            "z_score": round(lowess_result["z_score"], 4),
            "normalized_kiv2_coverage": round(cn_result["normalized_kiv2_coverage"], 4),
            "kiv2_repeat_copies": round(cn_result["kiv2_repeat_copies"], 4),
            "mad": round(mad, 4),
            "mad_flag": mad_flag,
        }

    except Exception as e:
        log(console, f"Error processing {sample_id}: {e}", style="danger")
        return None


# In[3]: Background region helpers
def _get_background_regions(
    span_start: int,
    span_end: int,
    excl_start: int,
    excl_end: int,
    window_size: int,
    num_windows: int,
    chrom: str = "chr6",
) -> list[tuple[str, int, int]]:
    starts = np.linspace(span_start, span_end - window_size, int(num_windows * 1.5), dtype=int)
    regions = []
    for s in starts:
        e = int(s) + window_size
        if excl_start <= s <= excl_end:
            continue
        regions.append((chrom, int(s), e))
        if len(regions) == num_windows:
            break
    return regions


def _compute_gc_for_regions(
    regions: list[tuple[str, int, int]],
    ref_fasta: str,
) -> dict[tuple[str, int, int], float]:
    gc_map = {}
    try:
        with pysam.FastaFile(ref_fasta) as fa:
            for chrom, start, end in regions:
                try:
                    seq = fa.fetch(chrom, start, end).upper()
                    total = len(seq) - seq.count("N")
                    gc = (
                        ((seq.count("G") + seq.count("C")) / total * 100)
                        if total > 0
                        else float("nan")
                    )
                    gc_map[(chrom, start, end)] = gc
                except Exception:
                    gc_map[(chrom, start, end)] = float("nan")
    except Exception:
        pass
    return gc_map


# In[4]: Read counting
def _count_reads(
    cram_path: str,
    ref_fasta: str,
    chrom: str,
    start: int,
    end: int,
    min_mapq: int,
) -> int:
    try:
        with pysam.AlignmentFile(cram_path, "rc", reference_filename=ref_fasta) as bam:
            refs = set(bam.references)
            target = chrom
            if target not in refs:
                alt = chrom[3:] if chrom.startswith("chr") else f"chr{chrom}"
                target = alt if alt in refs else None
            if target is None:
                return 0
            return bam.count(
                target,
                start,
                end,
                read_callback=lambda r: (
                    r.mapping_quality >= min_mapq and not r.is_secondary and not r.is_supplementary
                ),
            )
    except Exception:
        return 0


def _rpk(read_count: int, region_length: int) -> float:
    return (read_count / region_length) * 1000 if region_length > 0 else 0.0


# In[5]: LOWESS GC correction
def _lowess_gc_correct(
    lpa_gc: float,
    lpa_rpk: float,
    bg_df: pd.DataFrame,
    frac: float,
) -> dict:
    gc_vals = bg_df["gc"].values
    rpk_vals = bg_df["rpk"].values

    smoothed = lowess(rpk_vals, gc_vals, frac=frac, return_sorted=True)
    smooth_gc = smoothed[:, 0]
    smooth_rpk = smoothed[:, 1]

    predicted_bg = np.interp(gc_vals, smooth_gc, smooth_rpk)
    residuals_bg = rpk_vals - predicted_bg
    residual_mean = float(np.nanmean(residuals_bg))
    residual_std = float(np.nanstd(residuals_bg, ddof=1))

    predicted_lpa = float(np.interp(lpa_gc, smooth_gc, smooth_rpk))
    corrected_lpa = lpa_rpk / predicted_lpa if predicted_lpa > 0 else float("nan")
    residual_lpa = lpa_rpk - predicted_lpa
    z_score = (residual_lpa - residual_mean) / residual_std if residual_std > 0 else float("nan")

    return {
        "predicted_rpk": predicted_lpa,
        "corrected_rpk": corrected_lpa,
        "bg_residual_mean": residual_mean,
        "bg_residual_std": residual_std,
        "z_score": z_score,
    }


# In[6]: Copy number scaling (DRAGEN method)
def _compute_copy_number(
    corrected_rpk: float,
    mean_bg_rpk: float,
    extract_length: int,
    kiv2_ref_length: int,
    kiv2_scale_factor: float,
) -> dict:
    abs_corrected = corrected_rpk * mean_bg_rpk
    length_ratio = extract_length / kiv2_ref_length
    normalized_kiv2_coverage = abs_corrected * length_ratio
    kiv2_repeat_copies = normalized_kiv2_coverage * kiv2_scale_factor
    return {
        "normalized_kiv2_coverage": normalized_kiv2_coverage,
        "kiv2_repeat_copies": kiv2_repeat_copies,
    }


# In[7]: Thread-safe file writing
def _write_result(
    output_file: Path,
    result: dict,
    write_lock: Lock,
    header_written: list[bool],
) -> None:
    with write_lock:
        with open(output_file, "a") as f:
            if not header_written[0]:
                f.write("\t".join(result.keys()) + "\n")
                header_written[0] = True
            f.write("\t".join(str(v) for v in result.values()) + "\n")
