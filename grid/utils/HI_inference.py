from __future__ import annotations

import math
from collections import defaultdict
from pathlib import Path

from .utils import open_maybe_gz, log, progress_bar


def _read_dip_cn_file(dip_cn_file):
    """Read diploid CN file with string IDs, skipping non-data rows (e.g. header)."""
    IDs = []
    IRRs = []
    IDtoInd = {}
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
                continue  # skip header row
            id_ = parts[0]
            IDtoInd[id_] = len(IRRs)
            IDs.append(id_)
            IRRs.append(irr)
    return IDs, IRRs, IDtoInd


def _load_ibs_neighbors(neighbors_file, IDtoInd, MAX_NBR):
    """
    Load IBS neighbors from computeIBSpbwt output.

    File format (header + tab/space-separated rows):
        ID  hap  nbrInd  cMlen  cMedge  IDnbr  hapNbr

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
                h_idx = 2 * i + hap - 1  # 1-indexed to 0-indexed
                if len(hap_nbrs[h_idx]) < MAX_NBR:
                    hap_nbrs[h_idx].append((2 * j + hap_nbr - 1, 1.0))

    return hap_nbrs


def _segment_distance(bp1, bp2, region_start, region_end):
    """Physical distance from IBD segment [bp1, bp2] to target region. 0 if overlapping."""
    if bp2 < region_start:
        return float(region_start - bp2)
    if bp1 > region_end:
        return float(bp1 - region_end)
    return 0.0


def _load_ibd_neighbors(
    ilash_file,
    IDtoInd,
    MAX_NBR,
    region_start,
    region_end,
    min_length=0.5,
    min_match=0.70,
    weighted=False,
    weight_scale=1_000_000,
):
    """
    Load IBD neighbors from iLASH output.

    File format (no header, tab-separated):
        FID1  HAP_ID1  FID2  HAP_ID2  CHR  BP1  BP2  SNP_BP1  SNP_BP2  LENGTH  MATCH

    HAP_ID format: {FID}_{hap} where hap is 0 or 1 (0-indexed).
    LENGTH is in cM. MATCH is the IBD match score.

    If weighted=True, each neighbor gets a Lorentzian distance weight:
        w = (weight_scale / (distance_bp + weight_scale)) * match

    Returns hap_nbrs: list of lists of (neighbor_hap_idx, weight).
    """
    N = len(IDtoInd)
    # Collect all segments first so we can sort by length and truncate to MAX_NBR
    raw = defaultdict(list)  # h_idx -> [(nbr_h_idx, weight, length_cM)]

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
            hap_id1 = parts[1]  # e.g. HG00096_0
            fid2 = parts[2]
            hap_id2 = parts[3]  # e.g. HG00097_1

            try:
                bp1 = int(parts[5])
                bp2 = int(parts[6])
                length = float(parts[9])
                match = float(parts[10])
            except (ValueError, IndexError):
                continue

            if length < min_length or match < min_match:
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

            if weighted:
                dist = _segment_distance(bp1, bp2, region_start, region_end)
                w = (weight_scale / (dist + weight_scale)) * match
            else:
                w = 1.0

            h_idx1 = 2 * i + hap1
            h_idx2 = 2 * j + hap2
            raw[h_idx1].append((h_idx2, w, length))
            raw[h_idx2].append((h_idx1, w, length))

    hap_nbrs = [[] for _ in range(2 * N)]
    for h_idx, segments in raw.items():
        segments.sort(key=lambda x: -x[2])  # sort by length desc
        hap_nbrs[h_idx] = [(nbr, w) for nbr, w, _ in segments[:MAX_NBR]]

    return hap_nbrs


def _run_phasing(IRRs, hap_nbrs, MIN_NBR, N_ITERS, console):
    """
    Iterative haplotype phasing shared by both IBS and IBD methods.

    hap_nbrs[h] is a list of (neighbor_hap_idx, weight).
    Returns (hap_IRRs, mean_IRRs).
    """
    N = len(IRRs)
    hap_IRRs = [float("nan")] * (2 * N)

    N_to_phase = 0
    tot_nbrs = 0
    mean_IRRs = 0.0

    for i in range(N):
        h0 = hap_nbrs[2 * i]
        h1 = hap_nbrs[2 * i + 1]
        if len(h0) >= MIN_NBR and len(h1) >= MIN_NBR:
            hap_IRRs[2 * i] = IRRs[i] / 2
            hap_IRRs[2 * i + 1] = IRRs[i] / 2
            N_to_phase += 1
            tot_nbrs += len(h0) + len(h1)
            mean_IRRs += IRRs[i]

    log(console, f"Phasing {N_to_phase} samples with >={MIN_NBR} neighbors for both haps")
    if N_to_phase > 0:
        mean_IRRs /= N_to_phase
        log(console, f"Avg neighbors per hap: {tot_nbrs / (2.0 * N_to_phase):.2f}")

    for _ in range(N_ITERS):
        for i in range(N):
            if math.isnan(hap_IRRs[2 * i]):
                continue  # skip samples without enough neighbors

            hap_nbr_wsum = [1e-9, 1e-9]
            hap_nbr_wval = [0.0, 0.0]

            for h in range(2):
                for nbr, w in hap_nbrs[2 * i + h]:
                    val = hap_IRRs[nbr]
                    if not math.isnan(val):
                        hap_nbr_wsum[h] += w
                        hap_nbr_wval[h] += w * val

            nbr_mean_0 = hap_nbr_wval[0] / hap_nbr_wsum[0]
            nbr_mean_1 = hap_nbr_wval[1] / hap_nbr_wsum[1]
            denom = nbr_mean_0 + nbr_mean_1
            if denom > 0:
                hap_IRRs[2 * i] = IRRs[i] * nbr_mean_0 / denom
                hap_IRRs[2 * i + 1] = IRRs[i] * nbr_mean_1 / denom

    return hap_IRRs, mean_IRRs


def _compute_imp(i, hap_IRRs, hap_nbrs, mean_IRRs):
    """Compute imputed haplotype neighbor means for sample i at the final iteration."""
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

    # Fall back to global mean if no phased neighbors found
    if hap_nbr_wsum[0] <= 1e-9:
        imp0 = mean_IRRs / 2
    if hap_nbr_wsum[1] <= 1e-9:
        imp1 = mean_IRRs / 2

    return imp0, imp1


def hi_inference(config, console):
    """
    Impute haplotype-specific IRRs from diploid IRRs and haplotype neighbors.

    Supports two methods controlled by ``method`` in ``compute_haploid_genotypes``:

    * ``"ibs"`` — uses computeIBSpbwt output (``ibs_output`` key required)
    * ``"ibd"`` — uses iLASH output (``ibd_output`` key required).
      Set ``weighted: True`` to apply Lorentzian distance weighting.

    Output columns: ID  IRRs  hap1phased  hap2phased  hap1imp  hap2imp
    """
    try:
        hi_cfg = config.get("compute_haploid_genotypes", {})
        output_file_prefix = hi_cfg.get("output_file_prefix", "haploid_genotypes")
        output_file_type = config.get("output_file_type", "tsv")
        output_dir = config.get("output_dir", ".")
        output_file = Path(f"{output_dir}/{output_file_prefix}.{output_file_type}")

        dip_cn_file_prefix = config["compute_diploid_genotypes"].get("output_file_prefix")
        dip_cn_file = Path(f"{output_dir}/{dip_cn_file_prefix}.{output_file_type}")

        method = hi_cfg.get("method", "ibs").lower()
        MIN_NBR = hi_cfg.get("min_neighbors", 1)
        MAX_NBR = hi_cfg.get("max_neighbors", 10)
        N_ITERS = hi_cfg.get("n_iters", 100)
    except Exception as e:
        log(console, f"Config error: {e}", style="danger")
        return

    IDs, IRRs, IDtoInd = _read_dip_cn_file(dip_cn_file)
    N = len(IRRs)
    log(console, f"Read diploid IRR data for {N} samples", style="success")

    if method == "ibs":
        ibs_output = hi_cfg.get("ibs_output")
        if not ibs_output:
            log(console, "Config error: ibs_output required for method='ibs'", style="danger")
            return
        log(console, f"Loading IBS neighbors from {ibs_output}")
        hap_nbrs = _load_ibs_neighbors(ibs_output, IDtoInd, MAX_NBR)

    elif method == "ibd":
        ibd_output = hi_cfg.get("ibd_output")
        if not ibd_output:
            log(console, "Config error: ibd_output required for method='ibd'", style="danger")
            return
        region_start = config.get("start_bp")
        region_end = config.get("end_bp")
        weighted = hi_cfg.get("weighted", False)
        weight_scale = hi_cfg.get("weight_scale", 1_000_000)
        min_length = hi_cfg.get("min_length", 0.5)
        min_match = hi_cfg.get("min_match", 0.70)
        log(console, f"Loading IBD neighbors from {ibd_output} (weighted={weighted})")
        hap_nbrs = _load_ibd_neighbors(
            ibd_output,
            IDtoInd,
            MAX_NBR,
            region_start,
            region_end,
            min_length=min_length,
            min_match=min_match,
            weighted=weighted,
            weight_scale=weight_scale,
        )

    else:
        log(
            console,
            f"Config error: unknown method '{method}', must be 'ibs' or 'ibd'",
            style="danger",
        )
        return

    hap_IRRs, mean_IRRs = _run_phasing(IRRs, hap_nbrs, MIN_NBR, N_ITERS, console)

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
