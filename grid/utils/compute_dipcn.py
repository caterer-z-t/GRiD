# In[0]: Imports
import pandas as pd
import gzip
from pathlib import Path

from .utils import (
    log,
    progress_bar
)


# In[1]: Main function
def compute_diploid_genotypes(
        config,
        console
) -> None:
    """
    Replicate the awk normalization exactly:

        norm_reads = reads[sample] / sample_scale
                     / mean( reads[neighbor] / neighbor_scale )

    where the mean is over the top n_nbr neighbors.

    Args:
        config:  configuration dictionary
        console: Rich console for logging
    """
    try:
        output_file_prefix = config.get('compute_diploid_genotypes', {}).get('output_file_prefix', None)
        output_file_type   = config.get('output_file_type', 'tsv')
        output_dir         = config.get('output_dir', '.')
        output_file        = Path(f"{output_dir}/{output_file_prefix}.{output_file_type}")

        n_nbr = config.get('compute_diploid_genotypes', {}).get('n_nbr', 300)

        read_counts_file_prefix = config['count_reads'].get('output_file_prefix', None)
        read_counts_file = Path(f"{output_dir}/{read_counts_file_prefix}.{output_file_type}")

        zMax                  = config['mosdepth']['neighbors'].get('zmax', 2.0)
        neighbors_file_prefix = config['mosdepth']['neighbors'].get('output_file_prefix', None)
        neighbors_file        = Path(f"{output_dir}/{neighbors_file_prefix}.zMax{zMax:.1f}.{output_file_type}.gz")
    except Exception as e:
        log(console, f"Config error: {e}", style="danger")
        return

    # --- Load read counts ---
    reads_raw = pd.read_csv(read_counts_file, sep="\t", header=0, names=["Sample", "Reads"])
    reads_raw["Reads"] = pd.to_numeric(reads_raw["Reads"], errors="coerce")
    reads_raw.dropna(subset=["Reads"], inplace=True)
    reads: dict[str, float] = reads_raw.set_index("Sample")["Reads"].to_dict()

    # --- Load neighbors ---
    neighbors, sample_scales = load_neighbors(neighbors_file)

    # --- Compute normalized dipCN ---
    dipcn_list: list[tuple[str, float]] = []
    missing_ids: set[str] = set()

    with progress_bar(console, total=len(neighbors), description="Computing dipCN...") as (progress, task):
        for sample_id, nbr_list in neighbors.items():
            sample_scale = sample_scales.get(sample_id)
            if sample_scale is None or sample_id not in reads:
                progress.advance(task)
                continue

            sample_reads = reads[sample_id]

            total = 0.0
            count = 0
            for nbr_id, nbr_scale in nbr_list:
                if count >= n_nbr:
                    break
                if nbr_id not in reads:
                    missing_ids.add(nbr_id)
                    continue
                total += reads[nbr_id] / nbr_scale
                count += 1

            if count == 0:
                progress.advance(task)
                continue

            # Matches awk: reads[$1] / $2 / (sum / num)
            norm_reads = (sample_reads / sample_scale) / (total / count)
            dipcn_list.append((sample_id, norm_reads))
            progress.advance(task)

    if missing_ids:
        log(console,
            f"Warning: {len(missing_ids)} neighbor IDs not found in read counts "
            f"(showing up to 5: {list(missing_ids)[:5]})",
            style="warning")

    # --- Save ---
    dipcn_df = pd.DataFrame(dipcn_list, columns=["Sample", "Norm_Reads"])
    dipcn_df.to_csv(output_file, sep="\t", index=False)
    log(console, f"Saved {len(dipcn_df)} samples → {output_file}", style="success")


# In[2]: Load neighbors
def load_neighbors(neighbors_file: Path) -> tuple[dict[str, list[tuple[str, float]]], dict[str, float]]:
    """
    Parse the gzipped neighbors file.

    Each line format:
        sample_id  sample_scale  nbr1_id  nbr1_scale  nbr1_norm_dist  nbr2_id ...

    norm_dist is present but unused — only neighbor_id and neighbor_scale are needed.

    Returns:
        neighbors     : {sample_id: [(neighbor_id, neighbor_scale), ...]}
        sample_scales : {sample_id: sample_scale}
    """
    neighbors: dict[str, list[tuple[str, float]]] = {}
    sample_scales: dict[str, float] = {}

    with gzip.open(neighbors_file, "rt") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 2:
                continue

            sample_id = parts[0]
            try:
                sample_scale = float(parts[1])
            except ValueError:
                continue

            sample_scales[sample_id] = sample_scale

            nbr_list: list[tuple[str, float]] = []
            i = 2
            while i + 2 <= len(parts):
                nbr_id = parts[i]
                try:
                    nbr_scale = float(parts[i + 1])
                    # parts[i + 2] is norm_dist — unused
                except ValueError:
                    i += 3
                    continue
                nbr_list.append((nbr_id, nbr_scale))
                i += 3

            neighbors[sample_id] = nbr_list

    return neighbors, sample_scales