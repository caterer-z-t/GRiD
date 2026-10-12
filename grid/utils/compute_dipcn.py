# In[0]: Imports
import gzip
from pathlib import Path

from .utils import log, progress_bar


# In[1]: Main function
def compute_diploid_genotypes(config, console) -> None:
    """
    Compute diploid copy number estimates.

    Two modes depending on config:

    Mode 1 — GC normalize enabled (gc_normalize: run: True):
        Writes corrected_rpk_ratio directly from the GC normalization output
        to the diploid genotypes file. No neighbor normalization applied.
        This gives the best correlation with Bionano optical mapping (r=0.986).

    Mode 2 — GC normalize disabled (gc_normalize: run: False):
        Replicates Mukamel's awk normalization exactly:
            norm_reads = reads[sample] / sample_scale
                         / mean( reads[neighbor] / neighbor_scale )
        where the mean is over the top n_nbr neighbors.

    Args:
        config:  configuration dictionary
        console: Rich console for logging
    """
    try:
        output_file_prefix = config.get("compute_diploid_genotypes", {}).get(
            "output_file_prefix", "diploid_genotypes"
        )
        output_file_type = config.get("output_file_type", "tsv")
        output_dir = Path(config.get("output_dir", "."))
        output_file = output_dir / f"{output_file_prefix}.txt"
        n_nbr = config.get("compute_diploid_genotypes", {}).get("n_nbr", 300)

        read_counts_file_prefix = config["count_reads"].get("output_file_prefix", "read_counts")
        read_counts_file = output_dir / f"{read_counts_file_prefix}.{output_file_type}"

        zmax = config["mosdepth"]["neighbors"].get("zmax", 2.0)
        neighbors_file_prefix = config["mosdepth"]["neighbors"].get(
            "output_file_prefix", "ID_scale_neighbors"
        )
        neighbors_file = output_dir / f"{neighbors_file_prefix}.zMax{zmax:g}.txt.gz"

        gc_cfg = config.get("mosdepth", {}).get("gc_normalize", {})
        use_gc = gc_cfg.get("run", False)
        gc_prefix = gc_cfg.get("output_file_prefix", "gc_normalized")
        gc_file = output_dir / f"{gc_prefix}.{output_file_type}"

    except Exception as e:
        log(console, f"Config error: {e}", style="danger")
        return

    output_dir.mkdir(parents=True, exist_ok=True)

    # =========================================================================
    # Mode 1: GC normalize enabled — write corrected_rpk_ratio directly
    # =========================================================================
    if use_gc:
        if not gc_file.exists():
            log(console, f"GC normalized file not found: {gc_file}", style="danger")
            return

        log(
            console,
            "GC normalization enabled — writing corrected_rpk_ratio directly as diploid CN "
            "(no neighbor normalization)",
            style="info",
        )

        n_ok = 0
        n_na = 0
        with open(output_file, "w") as out_f:
            with open(gc_file, "r") as gc_f:
                header = gc_f.readline().strip().split("\t")
                try:
                    id_col = header.index("subject_id")
                    val_col = header.index("corrected_rpk_ratio")
                except ValueError:
                    log(
                        console,
                        "GC file missing 'subject_id' or 'corrected_rpk_ratio' columns.",
                        style="danger",
                    )
                    return

                for line in gc_f:
                    parts = line.strip().split("\t")
                    if len(parts) <= max(id_col, val_col):
                        continue
                    sample_id = parts[id_col]
                    try:
                        value = float(parts[val_col])
                        out_f.write(f"{sample_id} {value}\n")
                        n_ok += 1
                    except ValueError:
                        out_f.write(f"{sample_id} NA\n")
                        n_na += 1

        if n_na > 0:
            log(console, f"Warning: {n_na} samples written as NA", style="warning")
        log(
            console,
            f"Saved {n_ok} samples ({n_na} NA) → {output_file}",
            style="success",
        )
        return

    # =========================================================================
    # Mode 2: Raw read counts + neighbor normalization (Mukamel method)
    # =========================================================================
    if not read_counts_file.exists():
        log(console, f"Read counts file not found: {read_counts_file}", style="danger")
        return

    if not neighbors_file.exists():
        log(console, f"Neighbors file not found: {neighbors_file}", style="danger")
        return

    # Load raw read counts — no header, col1=ID col2=count
    reads: dict[str, float] = {}
    with open(read_counts_file, "r") as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 2:
                continue
            try:
                reads[parts[0]] = float(parts[1])
            except ValueError:
                continue
    log(console, f"Loaded raw read counts for {len(reads)} samples", style="info")

    # Mukamel neighbor normalization
    missing_ids: set[str] = set()
    n_na = 0
    n_ok = 0

    with open(output_file, "w") as out_f:
        with progress_bar(console, total=None, description="Computing dipCN...") as (
            progress,
            task,
        ):
            with gzip.open(neighbors_file, "rt") as nbr_f:
                for line in nbr_f:
                    parts = line.strip().split("\t")
                    if len(parts) < 2:
                        continue

                    sample_id = parts[0]
                    try:
                        sample_scale = float(parts[1])
                    except ValueError:
                        continue

                    total = 0.0
                    count = 0
                    i = 2
                    while i + 2 <= len(parts) and count < n_nbr:
                        nbr_id = parts[i]
                        try:
                            nbr_scale = float(parts[i + 1])
                        except ValueError:
                            i += 3
                            continue
                        if nbr_id in reads and nbr_scale > 0:
                            total += reads[nbr_id] / nbr_scale
                            count += 1
                        else:
                            missing_ids.add(nbr_id)
                        i += 3

                    if count > 0 and sample_scale > 0 and sample_id in reads:
                        norm_reads = (reads[sample_id] / sample_scale) / (total / count)
                        out_f.write(f"{sample_id} {norm_reads}\n")
                        n_ok += 1
                    else:
                        out_f.write(f"{sample_id} NA\n")
                        n_na += 1

                    progress.advance(task)

    if missing_ids:
        log(
            console,
            f"Warning: {len(missing_ids)} neighbor IDs not found in read counts "
            f"(showing up to 5: {list(missing_ids)[:5]})",
            style="warning",
        )
    if n_na > 0:
        log(console, f"Warning: {n_na} samples written as NA", style="warning")
    log(
        console,
        f"Saved {n_ok} samples ({n_na} NA) → {output_file}",
        style="success",
    )


# In[2]: Load neighbors (kept for potential reuse elsewhere)
def load_neighbors(
    neighbors_file: Path,
) -> tuple[dict[str, list[tuple[str, float]]], dict[str, float]]:
    """
    Parse the gzipped neighbors file.

    Each line format:
        sample_id  sample_scale  nbr1_id  nbr1_scale  nbr1_norm_dist  nbr2_id ...

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
                except ValueError:
                    i += 3
                    continue
                nbr_list.append((nbr_id, nbr_scale))
                i += 3

            neighbors[sample_id] = nbr_list

    return neighbors, sample_scales
