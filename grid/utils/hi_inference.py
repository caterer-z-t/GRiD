"""
hi_inference.py
===============
Haplotype-specific IRR inference for GRID.

Supports five phasing algorithms, all selectable via the config file key
``compute_haploid_genotypes.phasing_method``:

    iterative    — Hujoel C++ port: iteratively partition diploid IRR by
                   weighted neighbor means. The original default.

    graph_prop   — Label-spreading: blend own diploid/2 with neighbor mean
                   at each iteration via a tunable ``alpha``.
                   Best empirical performance (Pearson r ~0.88 vs. Bionano).

    oneshot      — Single-pass, no iteration. Fast; useful as a baseline.

    rank_based   — Like iterative but operates on rank-normalised hap IRRs,
                   reducing sensitivity to copy-number outliers.

    oneshot_rank — Single-pass rank-normalised variant of oneshot.

Two neighbor sources are supported (``method`` key):
    ibs  — computeIBSpbwt output
    ibd  — iLASH output

Example config snippet
----------------------
compute_haploid_genotypes:
  method: ibd                  # ibs | ibd
  phasing_method: graph_prop   # iterative | graph_prop | oneshot | rank_based | oneshot_rank
  ibd_output: path/to/output.match
  min_neighbors: 1
  max_neighbors: 20
  n_iters: 20                  # ignored by oneshot / oneshot_rank
  alpha: 0.7                   # graph_prop only
  agg_fn: mean                 # iterative | oneshot | rank_based | oneshot_rank
                               #   mean | median | weighted_cM | weighted_score |
                               #   weighted_cM_score | trimmed_mean
  # IBD-specific filters
  min_length: 2.0              # minimum segment length in cM
  min_match: 0.95              # minimum iLASH match score
  weighted: false              # Lorentzian distance weighting on neighbor load
  weight_scale: 1000000
"""

from __future__ import annotations

import math
import statistics
from collections import defaultdict
from pathlib import Path

from .utils import open_maybe_gz, log, progress_bar

def _read_dip_cn_file(dip_cn_file):
    """Read diploid CN file with string IDs, skipping non-data rows."""
    IDs, IRRs, IDtoInd = [], [], {}
    with open_maybe_gz(dip_cn_file) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            try:
                irr = float(parts[1])
            except ValueError:
                continue
            id_ = parts[0]
            IDtoInd[id_] = len(IRRs)
            IDs.append(id_)
            IRRs.append(irr)
    return IDs, IRRs, IDtoInd


def _load_ibs_neighbors(neighbors_file, IDtoInd, MAX_NBR):
    """
    Load IBS neighbors from computeIBSpbwt output.
    Haplotypes are 1-indexed (1 or 2).
    Returns hap_nbrs: list of lists of (neighbor_hap_idx, weight=1.0).
    """
    N = len(IDtoInd)
    hap_nbrs = [[] for _ in range(2 * N)]

    with open_maybe_gz(neighbors_file) as f:
        next(f)  # skip header
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) < 7:
                continue
            id_ = parts[0]
            try:
                hap = int(parts[1])
                hap_nbr = int(parts[6])
            except ValueError:
                continue
            id_nbr = parts[5]

            if hap not in (1, 2) or hap_nbr not in (1, 2):
                continue

            i = IDtoInd.get(id_)
            j = IDtoInd.get(id_nbr)
            if i is not None and j is not None:
                h_idx = 2 * i + hap - 1
                if len(hap_nbrs[h_idx]) < MAX_NBR:
                    hap_nbrs[h_idx].append((2 * j + hap_nbr - 1, 1.0))

    return hap_nbrs


def _adaptive_match_threshold(length_cM, min_length, min_match, match_at_min_length):
    """
    Length-dependent match score threshold.

    Short segments near min_length require a high match score
    (match_at_min_length).  Long segments asymptote down to min_match.

    Uses a hyperbolic curve fitted to two anchor points:
        (min_length,  match_at_min_length)
        (inf,         min_match)

    Derivation:
        threshold(L) = min_match + (match_at_min_length - min_match)
                       * min_length / L

    At L = min_length  →  match_at_min_length  (strict)
    At L → ∞           →  min_match            (lenient floor)

    If match_at_min_length == min_match the function is flat (disabled).
    """
    if match_at_min_length <= min_match:
        return min_match  # flat — adaptive mode disabled
    return min_match + (match_at_min_length - min_match) * min_length / length_cM


def _load_ibd_neighbors(
    ilash_file,
    IDtoInd,
    MAX_NBR,
    min_length=0.5,
    min_match=0.70,
    match_at_min_length=None,  # if None, falls back to flat min_match threshold
    region_start=None,
    region_end=None,
    weighted=False,
    weight_scale=1_000_000,
):
    """
    Load IBD neighbors from iLASH output using a two-tier priority system.

    Tier 1 (gold): segment overlaps the VNTR region [region_start, region_end]
                   AND passes the quality gate.
    Tier 2 (fill): passes quality gate but does not overlap the VNTR region.
                   Used to fill remaining neighbor slots after Tier 1.

    If region_start/region_end are None, all segments are Tier 1.

    Quality gate — length-adaptive match threshold:
        If match_at_min_length is set (and > min_match), the required match
        score scales with segment length:

            threshold(L) = min_match
                           + (match_at_min_length - min_match) * min_length / L

        Short segments (L ≈ min_length) must clear match_at_min_length.
        Long segments asymptote to min_match.
        If match_at_min_length is None or == min_match, a flat threshold
        (min_match) is used for all lengths.

    Within each tier, segments are ranked by length * match descending.
    HAP_ID format: {FID}_{hap} where hap is 0 or 1 (0-indexed).
    Returns hap_nbrs: list of lists of (neighbor_hap_idx, weight).
    """
    N = len(IDtoInd)
    use_tiers = region_start is not None and region_end is not None

    # raw[h_idx] = {"t1": [...], "t2": [...]}
    # each entry: (nbr_h_idx, weight, length_cM, match)
    raw = defaultdict(lambda: {"t1": [], "t2": []})

    with open_maybe_gz(ilash_file) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 11:
                parts = line.split()
            if len(parts) < 11:
                continue

            fid1 = parts[0]
            hap_id1 = parts[1]
            fid2 = parts[2]
            hap_id2 = parts[3]

            try:
                bp1 = int(parts[5])
                bp2 = int(parts[6])
                length = float(parts[9])
                match = float(parts[10])
            except (ValueError, IndexError):
                continue

            # quality gate — length-adaptive match threshold
            if length < min_length:
                continue
            effective_match = (
                _adaptive_match_threshold(length, min_length, min_match, match_at_min_length)
                if match_at_min_length is not None
                else min_match
            )
            if match < effective_match:
                continue

            try:
                hap1 = int(hap_id1.rsplit("_", 1)[-1])
                hap2 = int(hap_id2.rsplit("_", 1)[-1])
            except ValueError:
                continue

            if hap1 not in (0, 1) or hap2 not in (0, 1):
                continue

            i = IDtoInd.get(fid1)
            j = IDtoInd.get(fid2)
            if i is None or j is None:
                continue

            # assign tier
            if use_tiers:
                overlaps_region = (bp1 <= region_end) and (bp2 >= region_start)
                tier = "t1" if overlaps_region else "t2"
            else:
                tier = "t1"

            if weighted:
                w = (match**3) * length
            else:
                w = 1.0

            h_idx1 = 2 * i + hap1
            h_idx2 = 2 * j + hap2
            entry1 = (h_idx2, w, length, match)
            entry2 = (h_idx1, w, length, match)
            raw[h_idx1][tier].append(entry1)
            raw[h_idx2][tier].append(entry2)

    hap_nbrs = [[] for _ in range(2 * N)]
    sort_key = lambda x: -(x[2] * x[3])  # length * match descending

    t1_counts, t2_counts = [], []

    for h_idx in range(2 * N):
        buckets = raw[h_idx]
        t1_segs = sorted(buckets["t1"], key=sort_key)
        t2_segs = sorted(buckets["t2"], key=sort_key)

        # fill slots: Tier 1 first, then Tier 2 for remaining capacity
        chosen = t1_segs[:MAX_NBR]
        remaining = MAX_NBR - len(chosen)
        if remaining > 0:
            chosen += t2_segs[:remaining]

        hap_nbrs[h_idx] = [(nbr, w) for nbr, w, _, _ in chosen]
        t1_counts.append(len(t1_segs[:MAX_NBR]))
        t2_counts.append(len(t2_segs[:remaining]) if remaining > 0 else 0)

    if use_tiers:
        haps_with_t1 = sum(1 for c in t1_counts if c > 0)
        haps_with_any = sum(1 for h in hap_nbrs if len(h) > 0)
        avg_t1 = sum(t1_counts) / len(t1_counts) if t1_counts else 0
        avg_t2 = sum(t2_counts) / len(t2_counts) if t2_counts else 0
        return hap_nbrs, haps_with_t1, haps_with_any, avg_t1, avg_t2

    return hap_nbrs, None, None, None, None

def _collect_neighbor_vals(h_idx, hap_IRRs, hap_nbrs):
    """
    Return (vals, weights) for all non-NaN neighbors of haplotype h_idx.
    vals and weights are plain Python lists.
    """
    vals, weights = [], []
    for nbr, w in hap_nbrs[h_idx]:
        v = hap_IRRs[nbr]
        if not math.isnan(v):
            vals.append(v)
            weights.append(w)
    return vals, weights


def _aggregate(vals, weights, agg_fn, fallback):
    """
    Aggregate a list of neighbor values using the chosen function.

    agg_fn options
    --------------
    mean               — simple arithmetic mean (ignores weights)
    median             — median (ignores weights)
    weighted_cM        — weights proportional to segment cM  (stored in weights)
    weighted_score     — weights proportional to match score (stored in weights)
    weighted_cM_score  — combined cM × score weight         (stored in weights)
    trimmed_mean       — 10 % trimmed mean
    """
    if not vals:
        return fallback

    if agg_fn == "mean":
        return sum(vals) / len(vals)

    elif agg_fn == "median":
        return statistics.median(vals)

    elif agg_fn in ("weighted_cM", "weighted_score", "weighted_cM_score"):
        # For IBS neighbors weights are all 1.0, so these collapse to mean.
        # For IBD neighbors the weight field already encodes the chosen quantity
        # (set during neighbor loading). We just do a weighted average here.
        total_w = sum(weights)
        if total_w < 1e-12:
            return sum(vals) / len(vals)
        return sum(v * w for v, w in zip(vals, weights)) / total_w

    elif agg_fn == "trimmed_mean":
        n = len(vals)
        if n < 4:
            return sum(vals) / n
        cut = max(1, int(0.1 * n))
        trimmed = sorted(vals)[cut:-cut]
        return sum(trimmed) / len(trimmed)

    else:
        raise ValueError(
            f"Unknown agg_fn '{agg_fn}'. Choose from: mean, median, "
            "weighted_cM, weighted_score, weighted_cM_score, trimmed_mean"
        )

def _init_candidates(IRRs, hap_nbrs, MIN_NBR):
    """
    Identify samples with >= MIN_NBR neighbors on BOTH haplotypes.
    Returns (candidate_set, hap_IRRs initialised to diploid/2, mean_IRRs).
    """
    N = len(IRRs)
    hap_IRRs = [float("nan")] * (2 * N)
    candidates = set()
    mean_IRRs = 0.0

    for i in range(N):
        if len(hap_nbrs[2 * i]) >= MIN_NBR and len(hap_nbrs[2 * i + 1]) >= MIN_NBR:
            hap_IRRs[2 * i] = IRRs[i] / 2.0
            hap_IRRs[2 * i + 1] = IRRs[i] / 2.0
            candidates.add(i)
            mean_IRRs += IRRs[i]

    if candidates:
        mean_IRRs /= len(candidates)

    return candidates, hap_IRRs, mean_IRRs


# ── Algorithm 1: iterative (Hujoel C++ port) ─────────────────────────────────


def _phase_iterative(IRRs, hap_nbrs, MIN_NBR, N_ITERS, agg_fn, console):
    """
    Original Hujoel algorithm:
      initialise hap = diploid/2
      each iter: partition diploid by ratio of neighbor aggregate means
    """
    candidates, hap_IRRs, mean_IRRs = _init_candidates(IRRs, hap_nbrs, MIN_NBR)
    N = len(IRRs)

    log(
        console,
        f"[iterative] phasing {len(candidates)} samples, " f"agg_fn={agg_fn}, n_iters={N_ITERS}",
    )

    for _ in range(N_ITERS):
        new_IRRs = list(hap_IRRs)
        for i in candidates:
            means = []
            for h in range(2):
                vals, wts = _collect_neighbor_vals(2 * i + h, hap_IRRs, hap_nbrs)
                means.append(_aggregate(vals, wts, agg_fn, mean_IRRs / 2.0))
            denom = means[0] + means[1]
            if denom < 1e-9:
                denom = 1e-9
            new_IRRs[2 * i] = IRRs[i] * means[0] / denom
            new_IRRs[2 * i + 1] = IRRs[i] * means[1] / denom
        hap_IRRs = new_IRRs

    return hap_IRRs, mean_IRRs


# ── Algorithm 2: graph propagation (label-spreading) ─────────────────────────


def _phase_graph_prop(IRRs, hap_nbrs, MIN_NBR, N_ITERS, alpha, console):
    """
    Label-spreading:
      each iter: mixed_h = alpha*(diploid/2) + (1-alpha)*neighbor_mean
      then partition diploid by ratio of mixed values.

    alpha=0.7 (trust own signal 70 %) gave best empirical results.
    """
    candidates, hap_IRRs, mean_IRRs = _init_candidates(IRRs, hap_nbrs, MIN_NBR)
    N = len(IRRs)

    log(
        console,
        f"[graph_prop] phasing {len(candidates)} samples, " f"alpha={alpha}, n_iters={N_ITERS}",
    )

    for _ in range(N_ITERS):
        new_IRRs = list(hap_IRRs)
        for i in candidates:
            own_half = IRRs[i] / 2.0
            mixed = []
            for h in range(2):
                vals, wts = _collect_neighbor_vals(2 * i + h, hap_IRRs, hap_nbrs)
                nbr_mean = (sum(vals) / len(vals)) if vals else mean_IRRs / 2.0
                mixed.append(alpha * own_half + (1.0 - alpha) * nbr_mean)
            denom = mixed[0] + mixed[1]
            if denom < 1e-9:
                denom = 1e-9
            new_IRRs[2 * i] = IRRs[i] * mixed[0] / denom
            new_IRRs[2 * i + 1] = IRRs[i] * mixed[1] / denom
        hap_IRRs = new_IRRs

    return hap_IRRs, mean_IRRs


# ── Algorithm 3: oneshot ─────────────────────────────────────────────────────


def _phase_oneshot(IRRs, hap_nbrs, MIN_NBR, agg_fn, console):
    """
    Single-pass, non-iterative.
    Neighbours' values are always diploid/2 (the seed), so each sample's
    estimate is independent of all others.
    """
    N = len(IRRs)

    # seed: every hap starts at diploid/2
    seed_IRRs = []
    for irr in IRRs:
        seed_IRRs.extend([irr / 2.0, irr / 2.0])

    candidates = {
        i
        for i in range(N)
        if len(hap_nbrs[2 * i]) >= MIN_NBR and len(hap_nbrs[2 * i + 1]) >= MIN_NBR
    }

    mean_IRRs = sum(IRRs[i] for i in candidates) / len(candidates) if candidates else 1.0

    hap_IRRs = [float("nan")] * (2 * N)
    log(console, f"[oneshot] phasing {len(candidates)} samples, agg_fn={agg_fn}")

    for i in candidates:
        means = []
        for h in range(2):
            vals, wts = _collect_neighbor_vals(2 * i + h, seed_IRRs, hap_nbrs)
            means.append(_aggregate(vals, wts, agg_fn, mean_IRRs / 2.0))
        denom = means[0] + means[1]
        if denom < 1e-9:
            denom = 1e-9
        hap_IRRs[2 * i] = IRRs[i] * means[0] / denom
        hap_IRRs[2 * i + 1] = IRRs[i] * means[1] / denom

    return hap_IRRs, mean_IRRs


# ── Algorithm 4: rank-based iterative ────────────────────────────────────────


def _phase_rank_based(IRRs, hap_nbrs, MIN_NBR, N_ITERS, agg_fn, console):
    """
    Iterative, but ranks hap IRRs to [0,1] before aggregating neighbours.
    Reduces influence of extreme copy-number outliers on the mean signal.
    """
    candidates, hap_IRRs, mean_IRRs = _init_candidates(IRRs, hap_nbrs, MIN_NBR)
    N = len(IRRs)

    log(
        console,
        f"[rank_based] phasing {len(candidates)} samples, " f"agg_fn={agg_fn}, n_iters={N_ITERS}",
    )

    phased_indices = [2 * i + h for i in candidates for h in range(2)]

    for _ in range(N_ITERS):
        # rank-normalise current hap IRRs for phased haplotypes only
        raw_vals = [hap_IRRs[idx] for idx in phased_indices]
        sorted_vals = sorted(raw_vals)
        n_phased = len(sorted_vals)
        rank_map = {
            idx: (sorted_vals.index(hap_IRRs[idx]) + 1) / n_phased for idx in phased_indices
        }

        # build a ranked version for neighbor lookup
        ranked_IRRs = list(hap_IRRs)
        for idx, rank in rank_map.items():
            ranked_IRRs[idx] = rank

        new_IRRs = list(hap_IRRs)
        for i in candidates:
            means = []
            for h in range(2):
                vals, wts = _collect_neighbor_vals(2 * i + h, ranked_IRRs, hap_nbrs)
                means.append(_aggregate(vals, wts, agg_fn, 0.5))
            denom = means[0] + means[1]
            if denom < 1e-9:
                denom = 1e-9
            new_IRRs[2 * i] = IRRs[i] * means[0] / denom
            new_IRRs[2 * i + 1] = IRRs[i] * means[1] / denom
        hap_IRRs = new_IRRs

    return hap_IRRs, mean_IRRs


# ── Algorithm 5: oneshot rank ─────────────────────────────────────────────────


def _phase_oneshot_rank(IRRs, hap_nbrs, MIN_NBR, agg_fn, console):
    """
    Single-pass rank-normalised variant.
    Seeds are rank-normalised diploid/2 values; no iteration.
    """
    N = len(IRRs)

    half_vals = [irr / 2.0 for irr in IRRs]
    sorted_half = sorted(half_vals)
    n_all = len(sorted_half)
    seed_IRRs = []
    for v in half_vals:
        rank = (sorted_half.index(v) + 1) / n_all
        seed_IRRs.extend([rank, rank])

    candidates = {
        i
        for i in range(N)
        if len(hap_nbrs[2 * i]) >= MIN_NBR and len(hap_nbrs[2 * i + 1]) >= MIN_NBR
    }

    mean_IRRs = sum(IRRs[i] for i in candidates) / len(candidates) if candidates else 1.0

    hap_IRRs = [float("nan")] * (2 * N)
    log(console, f"[oneshot_rank] phasing {len(candidates)} samples, agg_fn={agg_fn}")

    for i in candidates:
        means = []
        for h in range(2):
            vals, wts = _collect_neighbor_vals(2 * i + h, seed_IRRs, hap_nbrs)
            means.append(_aggregate(vals, wts, agg_fn, 0.5))
        denom = means[0] + means[1]
        if denom < 1e-9:
            denom = 1e-9
        hap_IRRs[2 * i] = IRRs[i] * means[0] / denom
        hap_IRRs[2 * i + 1] = IRRs[i] * means[1] / denom

    return hap_IRRs, mean_IRRs

def _compute_imp(i, hap_IRRs, hap_nbrs, mean_IRRs):
    """Compute imputed haplotype neighbour means for sample i."""
    hap_nbr_wsum = [1e-9, 1e-9]
    hap_nbr_wval = [0.0, 0.0]

    for h in range(2):
        for nbr, w in hap_nbrs[2 * i + h]:
            val = hap_IRRs[nbr]
            if not math.isnan(val):
                hap_nbr_wsum[h] += w
                hap_nbr_wval[h] += w * val

    imp0 = hap_nbr_wval[0] / hap_nbr_wsum[0]
    imp1 = hap_nbr_wval[1] / hap_nbr_wsum[1]

    if hap_nbr_wsum[0] <= 1e-9:
        imp0 = mean_IRRs / 2.0
    if hap_nbr_wsum[1] <= 1e-9:
        imp1 = mean_IRRs / 2.0

    return imp0, imp1

#: All valid values for phasing_method and agg_fn (for validation / help text)
PHASING_METHODS = ("iterative", "graph_prop", "oneshot", "rank_based", "oneshot_rank")
AGG_FNS = ("mean", "median", "weighted_cM", "weighted_score", "weighted_cM_score", "trimmed_mean")


def hi_inference(config, console):
    """
    Impute haplotype-specific IRRs from diploid IRRs and haplotype neighbours.

    Neighbour source  (``method`` key)
    -----------------------------------
    ibs   — computeIBSpbwt output (``ibs_output`` required)
    ibd   — iLASH output          (``ibd_output`` required)

    Phasing algorithm  (``phasing_method`` key, default: ``graph_prop``)
    ---------------------------------------------------------------------
    iterative     — Hujoel C++ port (original GRID default)
    graph_prop    — Label-spreading; tune with ``alpha`` (default 0.7)
    oneshot       — Single-pass, no iteration
    rank_based    — Iterative on rank-normalised values
    oneshot_rank  — Single-pass rank-normalised

    Aggregation function  (``agg_fn`` key, default: ``mean``)
    ----------------------------------------------------------
    Used by iterative, oneshot, rank_based, oneshot_rank.
    Ignored by graph_prop (always uses plain mean for neighbour signal).

        mean | median | weighted_cM | weighted_score | weighted_cM_score | trimmed_mean

    Output columns: ID  IRRs  hap1phased  hap2phased  hap1imp  hap2imp
    """
    # ── Parse config ──────────────────────────────────────────────────────────
    try:
        hi_cfg = config.get("compute_haploid_genotypes", {})

        output_file_prefix = hi_cfg.get("output_file_prefix", "haploid_genotypes")
        output_file_type = config.get("output_file_type", "tsv")
        output_dir = config.get("output_dir", ".")
        output_file = Path(f"{output_dir}/{output_file_prefix}.{output_file_type}")

        dip_cn_file_prefix = config["compute_diploid_genotypes"].get("output_file_prefix")
        # The diploid step always writes .txt regardless of output_file_type,
        # so check the configured extension first, then fall back through common
        # variants to avoid a silent "file not found" crash.
        _candidates = [output_file_type, "txt", "tsv", "txt.gz", "tsv.gz"]
        dip_cn_file = None
        for ext in _candidates:
            candidate = Path(f"{output_dir}/{dip_cn_file_prefix}.{ext}")
            if candidate.exists():
                dip_cn_file = candidate
                break
        if dip_cn_file is None:
            log(
                console,
                f"Diploid genotypes file not found under {output_dir}/ — "
                f"tried extensions: {_candidates}. "
                f"Ensure compute_diploid_genotypes ran successfully.",
                style="danger",
            )
            return

        method = hi_cfg.get("method", "ibs").lower()
        phasing_method = hi_cfg.get("phasing_method", "graph_prop").lower()
        MIN_NBR = hi_cfg.get("min_neighbors", 1)
        MAX_NBR = hi_cfg.get("max_neighbors", 20)
        N_ITERS = hi_cfg.get("n_iters", 20)
        alpha = hi_cfg.get("alpha", 0.7)
        agg_fn = hi_cfg.get("agg_fn", "mean").lower()

    except Exception as e:
        log(console, f"Config error: {e}", style="danger")
        return

    # Validate user choices
    if phasing_method not in PHASING_METHODS:
        log(
            console,
            f"Unknown phasing_method '{phasing_method}'. "
            f"Choose from: {', '.join(PHASING_METHODS)}",
            style="danger",
        )
        return
    if agg_fn not in AGG_FNS:
        log(
            console,
            f"Unknown agg_fn '{agg_fn}'. " f"Choose from: {', '.join(AGG_FNS)}",
            style="danger",
        )
        return

    # ── Load diploid data ─────────────────────────────────────────────────────
    IDs, IRRs, IDtoInd = _read_dip_cn_file(dip_cn_file)
    N = len(IRRs)
    log(console, f"Read diploid IRR data for {N} samples", style="success")

    # ── Load neighbour graph ──────────────────────────────────────────────────
    if method == "ibs":
        ibs_output = hi_cfg.get("ibs_output")
        if not ibs_output:
            log(console, "Config error: ibs_output required for method='ibs'", style="danger")
            return
        log(console, f"Loading IBS neighbours from {ibs_output}")
        hap_nbrs = _load_ibs_neighbors(ibs_output, IDtoInd, MAX_NBR)

    elif method == "ibd":
        ibd_output = hi_cfg.get("ibd_output")
        if not ibd_output:
            log(console, "Config error: ibd_output required for method='ibd'", style="danger")
            return
        weighted = hi_cfg.get("weighted", False)
        weight_scale = hi_cfg.get("weight_scale", 1_000_000)
        min_length = hi_cfg.get("min_length", 0.5)
        min_match = hi_cfg.get("min_match", 0.70)
        match_at_min_length = hi_cfg.get("match_at_min_length", None)
        region_start = config.get("start_bp", None)
        region_end = config.get("end_bp", None)

        use_tiers = region_start is not None and region_end is not None
        use_adaptive = match_at_min_length is not None and match_at_min_length > min_match

        if use_adaptive:
            log(
                console,
                f"Adaptive match threshold: "
                f"{match_at_min_length:.3f} at {min_length} cM → "
                f"{min_match:.3f} at long segments  "
                f"[threshold(L) = {min_match:.3f} + "
                f"{match_at_min_length - min_match:.3f} × {min_length}/L]",
            )
        log(
            console,
            f"Loading IBD neighbours from {ibd_output} "
            f"(weighted={weighted}, min_length={min_length}, "
            f"min_match={'adaptive' if use_adaptive else min_match}, "
            f"tiered={'yes' if use_tiers else 'no'})",
        )

        hap_nbrs, haps_t1, haps_any, avg_t1, avg_t2 = _load_ibd_neighbors(
            ibd_output,
            IDtoInd,
            MAX_NBR,
            min_length=min_length,
            min_match=min_match,
            match_at_min_length=match_at_min_length,
            region_start=region_start,
            region_end=region_end,
            weighted=weighted,
            weight_scale=weight_scale,
        )

        if use_tiers:
            log(
                console,
                f"Tier stats (per haplotype): "
                f"{haps_t1} haps have ≥1 Tier-1 (VNTR-overlapping) neighbour | "
                f"{haps_any} haps have ≥1 neighbour of any tier | "
                f"avg Tier-1 nbrs: {avg_t1:.2f}  avg Tier-2 nbrs: {avg_t2:.2f}",
            )

    else:
        log(console, f"Unknown method '{method}'. Must be 'ibs' or 'ibd'.", style="danger")
        return

    # ── Run chosen phasing algorithm ──────────────────────────────────────────
    log(
        console,
        f"Running phasing_method='{phasing_method}' "
        f"(agg_fn='{agg_fn}', alpha={alpha}, n_iters={N_ITERS})",
    )

    if phasing_method == "iterative":
        hap_IRRs, mean_IRRs = _phase_iterative(IRRs, hap_nbrs, MIN_NBR, N_ITERS, agg_fn, console)

    elif phasing_method == "graph_prop":
        hap_IRRs, mean_IRRs = _phase_graph_prop(IRRs, hap_nbrs, MIN_NBR, N_ITERS, alpha, console)

    elif phasing_method == "oneshot":
        hap_IRRs, mean_IRRs = _phase_oneshot(IRRs, hap_nbrs, MIN_NBR, agg_fn, console)

    elif phasing_method == "rank_based":
        hap_IRRs, mean_IRRs = _phase_rank_based(IRRs, hap_nbrs, MIN_NBR, N_ITERS, agg_fn, console)

    elif phasing_method == "oneshot_rank":
        hap_IRRs, mean_IRRs = _phase_oneshot_rank(IRRs, hap_nbrs, MIN_NBR, agg_fn, console)

    # ── Write output ──────────────────────────────────────────────────────────
    with open_maybe_gz(output_file, "wt") as fout:
        fout.write("ID\tIRRs\thap1phased\thap2phased\thap1imp\thap2imp\n")
        for i in range(N):
            imp0, imp1 = _compute_imp(i, hap_IRRs, hap_nbrs, mean_IRRs)
            fout.write(
                f"{IDs[i]}\t{IRRs[i]:.2f}\t"
                f"{hap_IRRs[2*i]:.2f}\t{hap_IRRs[2*i+1]:.2f}\t"
                f"{imp0:.2f}\t{imp1:.2f}\n"
            )

    log(console, f"Haploid genotypes written to {output_file}", style="success")
