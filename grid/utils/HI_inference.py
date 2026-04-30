# In[0]: Imports
import math
from pathlib import Path

from .utils import (
    open_maybe_gz,
    log, 
    progress_bar
)


def read_dip_cn_file(dip_cn_file, IDs, IRRs, IDtoInd):
    """
    Reads the diploid copy number file and populates the IRRs, IDs, and IDtoInd structures.

    Args:
        dip_cn_file: Path to input file containing sample-level IRRs (ID, IRR)
        IDs: list to populate with sample IDs
        IRRs: list to populate with sample-level IRRs
        IDtoInd: dict to populate mapping sample ID to index in IRRs and IDs lists

    Returns:        
        IRRs: list of sample-level IRRs
        IDs: list of sample IDs
        IDtoInd: dict mapping sample ID to index in IRRs and IDs lists
    """
    with open_maybe_gz(dip_cn_file) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            id_ = int(parts[0])
            irr = float(parts[1])
            IDtoInd[id_] = len(IRRs)
            IDs.append(id_)
            IRRs.append(irr)
    return IRRs, IDs, IDtoInd

def read_neighbors_file(neighbors_file, IDtoInd, hap_nbrs, MAX_NBR):
    """
    Reads the neighbors file and populates the hap_nbrs structure.

    Args:
        neighbors_file: Path to input file containing haplotype neighbor info (ID, hap, nbr_ID, ...)
        IDtoInd: dict mapping sample ID to index in IRRs and IDs lists
        hap_nbrs: list to populate with neighbor indices for each haplotype 
        MAX_NBR: maximum number of neighbors to store per haplotype
    Returns:
        hap_nbrs: list of neighbor indices for each haplotype
    """
    with open_maybe_gz(neighbors_file) as f:
        next(f)  # skip header
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            id_, hap, nbr_ind = int(parts[0]), int(parts[1]), int(parts[2])
            cM_len, cM_edge = float(parts[3]), float(parts[4])
            id_nbr, hap_nbr = int(parts[5]), int(parts[6])

            if id_ < 0 or id_nbr < 0:
                continue

            assert hap in (1, 2)

            i = IDtoInd.get(id_)
            j = IDtoInd.get(id_nbr)
            if i is not None and j is not None:
                h_idx = 2 * i + hap - 1
                if len(hap_nbrs[h_idx]) < MAX_NBR:
                    hap_nbrs[h_idx].append(2 * j + hap_nbr - 1)

    return hap_nbrs


def hi_inference(config, console):
    """
    Imputes haplotype-specific IRRs from sample-level IRRs and haplotype neighbors, using an iterative approach.

    Args:
        dip_cn_file: Path to input file containing sample-level IRRs (ID, IRR)
        neighbors_file: Path to input file containing haplotype neighbor info (ID, hap, nbr_ID, ...)
        output_file: Path to output file for imputed haplotype-specific IRRs
        console: Rich console for logging

    Output file format:
    ID  IRRs    hap1phased  hap2phased  hap1imp    hap2imp
     - ID: sample ID
     - IRRs: original sample-level IRR
     - hap1phased, hap2phased: phased haplotype-specific IRRs after iterative phasing
     - hap1imp, hap2imp: mean-imputed haplotype-specific IRRs based on neighbors (for samples with no phased neighbors)
    """
    try:
        output_file_prefix = config.get('HI_inference', {}).get('output_file_prefix', None)
        output_file_type   = config.get('output_file_type', 'tsv')
        output_dir         = config.get('output_dir', '.')
        output_file        = Path(f"{output_dir}/{output_file_prefix}.{output_file_type}")

        dip_cn_file_prefix = config['compute_diploid_genotypes'].get('output_file_prefix', None)
        dip_cn_file = Path(f"{output_dir}/{dip_cn_file_prefix}.{output_file_type}")

        zMax                  = config['mosdepth']['neighbors'].get('zmax', 2.0)
        neighbors_file_prefix = config['mosdepth']['neighbors'].get('output_file_prefix', None)
        neighbors_file        = Path(f"{output_dir}/{neighbors_file_prefix}.zMax{zMax:.1f}.{output_file_type}.gz")

        MIN_NBR = config.get('HI_inference', {}).get('min_neighbors', 1)
        MAX_NBR = config.get('HI_inference', {}).get('max_neighbors', 10)
        N_ITERS = config.get('HI_inference', {}).get('n_iters', 20)

    except Exception as e:
        log(console, f"Config error: {e}", style="danger")
        return

    IRRs, IDs, IDtoInd = read_dip_cn_file(dip_cn_file, [], [], {})

    N = len(IRRs)
    log(console, f"Read IRR data for {N} samples", style="success")

    hap_nbrs = [[] for _ in range(2 * N)]  # hap_nbrs[2*i], hap_nbrs[2*i+1] for sample i
    hap_IRRs = [float("nan")] * (2 * N)

    hap_nbrs = read_neighbors_file(neighbors_file, IDtoInd, hap_nbrs, MAX_NBR)

    # --- initialize hapIRRs to half of IRRs ---
    N_to_phase = 0
    tot_nbrs = 0
    mean_IRRs = 0.0

    for i in range(N):
        if len(hap_nbrs[2 * i]) >= MIN_NBR and len(hap_nbrs[2 * i + 1]) >= MIN_NBR:
            hap_IRRs[2 * i] = IRRs[i] / 2
            hap_IRRs[2 * i + 1] = IRRs[i] / 2
            N_to_phase += 1
            tot_nbrs += len(hap_nbrs[2 * i]) + len(hap_nbrs[2 * i + 1])
            mean_IRRs += IRRs[i]

    log(console, f"Read hap neighbors; phasing {N_to_phase} samples with >={MIN_NBR} neighbors for both haps")
    if N_to_phase > 0:
        log(console, f"Average number of neighbors to use per hap: {tot_nbrs / (2.0 * N_to_phase):.4f}")
        mean_IRRs /= N_to_phase

    with open_maybe_gz(output_file, "wt") as fout:
        fout.write("ID\tIRRs\thap1phased\thap2phased\thap1imp\thap2imp\n")

        for iter_ in range(N_ITERS):
            last_iter = (iter_ == N_ITERS - 1)

            for i in range(N):
                hap_nbr_sum = [1e-9, 1e-9]
                hap_nbr_num = [0, 0]
                hap_nbr_mean = [0.0, 0.0]

                for h in range(2):
                    for nbr in hap_nbrs[2 * i + h]:
                        val = hap_IRRs[nbr]
                        if not math.isnan(val):
                            hap_nbr_sum[h] += val
                            hap_nbr_num[h] += 1
                    if hap_nbr_num[h]:
                        hap_nbr_mean[h] = hap_nbr_sum[h] / hap_nbr_num[h]
                    else:
                        hap_nbr_mean[h] = mean_IRRs / 2  # mean-impute if no neighbors

                denom = hap_nbr_mean[0] + hap_nbr_mean[1]
                for h in range(2):
                    if not math.isnan(hap_IRRs[2 * i + h]):
                        hap_IRRs[2 * i + h] = IRRs[i] * hap_nbr_mean[h] / denom

                if last_iter:
                    fout.write(
                        f"{IDs[i]}\t{IRRs[i]:.2f}\t"
                        f"{hap_IRRs[2*i]:.2f}\t{hap_IRRs[2*i+1]:.2f}\t"
                        f"{hap_nbr_mean[0]:.2f}\t{hap_nbr_mean[1]:.2f}\n"
                    )