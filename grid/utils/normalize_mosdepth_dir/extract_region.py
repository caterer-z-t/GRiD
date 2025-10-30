# grid/utils/normalize_mosdepth_dir/extract_region.py
# In[1]: Imports
import gzip
from pathlib import Path
from collections import defaultdict

from .get_regions_bed import regions_bed_gz

# In[2]: Function to Extract Regions
def extract_regions(
    individuals: dict[str, Path],
    chromosome: str,
    start: int,
    end: int,
    min_depth: int = 20,
    max_depth: int = 100,
    excluded: dict[str, set[int]] = None,
    progress=None,
    task=None
) -> dict[str, list[tuple[int, int, float]]]:
    """
    Extract depth regions passing filters and not overlapping repeats.

    Args:
        individuals: dict mapping sample IDs to their mosdepth.global.dist file paths
        chromosome: Chromosome name to extract
        start: Start position of region
        end: End position of region
        min_depth: Minimum depth threshold
        max_depth: Maximum depth threshold
        excluded: dict mapping chromosome to set of 1kb bins to exclude (from repeat mask)
        progress: Optional Rich Progress object for updates
        task: Optional Rich Progress task ID for updates

    Returns:
        Dictionary mapping {sample_id: [(start, end, depth), ...]} for extracted regions
    """
    regions_to_extract = defaultdict(list)

    # get the integer from the chromosome string if it starts with 'chr'
    if chromosome.startswith("chr"):
        chromosome = chromosome[3:]

    for ind_id, global_dist_file in individuals.items():
        bed_gz = regions_bed_gz(global_dist_file)

        if not bed_gz.exists():
            if progress:
                progress.console.print(f"[yellow]Warning:[/yellow] Missing {bed_gz} for {ind_id}")
            if task:
                progress.update(task, advance=1)
            continue

        try:
            with gzip.open(bed_gz, "rt") as f:
                for line in f:
                    if not line.startswith(chromosome):
                        continue
                    fields = line.strip().split("\t")
                    if len(fields) < 4:
                        continue
                    chrom, r_start, r_end, depth = (
                        fields[0], int(fields[1]), int(fields[2]), float(fields[3])
                    )

                    # Skip if no depth or outside region
                    if not (depth > 0 and r_end >= start and r_start <= end):
                        continue

                    # Skip if overlaps repeat mask
                    if excluded:
                        region_kb = set(range(r_start // 1000, r_end // 1000 + 1))
                        if region_kb & excluded.get(chrom, set()):
                            continue

                    # Skip if outside depth range
                    if depth < min_depth or depth > max_depth:
                        continue

                    regions_to_extract[ind_id].append((r_start, r_end, depth))

        except OSError as e:
            # Handle truncated or uncompressed files gracefully
            if progress:
                progress.console.print(f"[red]Error reading {bed_gz} for {ind_id}: {e}[/red]")
            continue

        finally:
            if task:
                progress.update(task, advance=1)

    return regions_to_extract
