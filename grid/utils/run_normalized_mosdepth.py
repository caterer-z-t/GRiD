# grid/utils/run_normalize_mosdepth.py
# In[1]: Imports
from pathlib import Path
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed

from rich.console import Console
from rich.progress import Progress, SpinnerColumn, TextColumn, BarColumn, TaskProgressColumn

from .normalize_mosdepth_dir.get_individuals import get_individuals
from .normalize_mosdepth_dir.load_repeat_mask import load_repeat_mask
from .normalize_mosdepth_dir.process_one_individual import _process_one_individual
from .normalize_mosdepth_dir.build_depth_matrix import build_matrix_from_regions
from .normalize_mosdepth_dir.normalize_matrix import normalize_matrix
from .normalize_mosdepth_dir.select_high_variance_region import select_high_variance_regions
from .normalize_mosdepth_dir.write_normalized_output import write_normalized_output
from .helper_dir.create_region import create_region_string


# In[2]: Main Function to Run Normalize Mosdepth
def run_normalize_mosdepth(
    mosdepth_dir: str,
    output_file: str,
    repeat_mask: str,
    chrom: str,
    start: int,
    end: int,
    min_depth: int = 20,
    max_depth: int = 100,
    top_frac: float = 0.1,
    threads: int = 1,
    console: Console = None,
):
    """
    Normalize mosdepth coverage across all samples in a directory.
    
    Args:
        mosdepth_dir: Directory containing mosdepth output files
        output_file: Path to output normalized .tsv.gz file
        repeat_mask: Path to repeat mask file
        chrom: Chromosome name
        start: Start position
        end: End position
        min_depth: Minimum depth threshold for regions
        max_depth: Maximum depth threshold for regions
        top_frac: Fraction of high-variance regions to select
        threads: Number of threads for parallel processing

    Returns:
        None
    """
    if console is None:
        console = Console()

    output_path = Path(output_file).expanduser()
    output_path.parent.mkdir(parents=True, exist_ok=True)

    console.rule("[bold blue]Step 1: Collect Individuals")
    individuals = get_individuals(mosdepth_dir)
    if not individuals:
        console.print(f"[red]✗ No mosdepth files found in {mosdepth_dir}[/red]")
        sys.exit(1)
    console.print(f"[green]✓ Found {len(individuals)} samples[/green]")

    console.rule("[bold blue]Step 2: Load Repeat Mask")
    excluded = load_repeat_mask(repeat_mask)
    console.print(f"[green]✓ Loaded {sum(len(v) for v in excluded.values()):,} excluded bins[/green]")

    console.rule("[bold blue]Step 3: Extract Regions")
    region_str = create_region_string(None, chrom, start, end)
    console.print(f"[yellow]Extracting regions for {region_str} using {threads} thread(s)[/yellow]")

    regions_to_extract = {}
    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        TaskProgressColumn(),
        console=console,
    ) as progress:
        task = progress.add_task("[green]Extracting per-sample regions...", total=len(individuals))

        with ThreadPoolExecutor(max_workers=max(1, threads)) as executor:
            futures = {
                executor.submit(
                    _process_one_individual,
                    (ind_id, global_file, chrom, start, end, min_depth, max_depth, excluded),
                ): ind_id
                for ind_id, global_file in individuals.items()
            }

            for future in as_completed(futures):
                ind_id = futures[future]
                try:
                    ind_id, regions = future.result()
                    regions_to_extract[ind_id] = regions
                except Exception as e:
                    console.print(f"[red]Error processing {ind_id}: {e}[/red]")
                finally:
                    progress.update(task, advance=1)

    n_regions = sum(len(v) for v in regions_to_extract.values())
    console.print(f"[green]✓ Extracted {n_regions:,} total valid regions[/green]")

    console.rule("[bold blue]Step 4: Build Depth Matrix")
    individuals_order, regions_list, mat = build_matrix_from_regions(regions_to_extract)
    console.print(f"[green]✓ Matrix built: {mat.shape[0]} samples × {mat.shape[1]} regions[/green]")

    console.rule("[bold blue]Step 5: Normalize Matrix")
    normalized_mat, variance_ratios = normalize_matrix(mat)
    console.print(f"[green]✓ Normalization complete[/green]")

    console.rule("[bold blue]Step 6: Select High-Variance Regions")
    selected_indices = select_high_variance_regions(variance_ratios, top_frac)
    console.print(f"[green]✓ Selected top {len(selected_indices)} regions ({top_frac*100:.1f}%)[/green]")

    console.rule("[bold blue]Step 7: Write Output")
    write_normalized_output(
        normalized_mat,
        individuals_order,
        regions_list,
        selected_indices,
        output_path,
    )
    console.print(f"[bold blue]✓ Done. Output written to {output_path}[/bold blue]")
