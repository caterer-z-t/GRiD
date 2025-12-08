from pathlib import Path
from .subset_cram import subset_cram
from concurrent.futures import ThreadPoolExecutor, as_completed
from threading import Lock
from functools import partial
from rich.console import Console
from rich.progress import Progress, SpinnerColumn, TextColumn, BarColumn, TaskProgressColumn

def batch_subset_cram(cram_dir, region, output_dir, reference, console):
    """
    Subset all CRAM files in a directory for a specific genomic region.

    Args:
        cram_dir: Directory containing CRAM files.
        region: Region in 'chr:start-end' format (e.g., 'chr6:160000000-160100000').
        output_dir: Directory to save subset CRAM files.
        reference: Reference genome FASTA (required for CRAM output).
        console: Console object for logging.
        
    Returns:
        List of paths to the subset CRAM files.
    """

    cram_path = Path(cram_dir)
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    write_lock = Lock()

    process_func = partial(
        subset_cram,
        region=region,
        reference=reference
    )

    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        TaskProgressColumn(),
        console=console
    ) as progress:

        task = progress.add_task("[bold green] Subsetting CRAM files...", total=len(cram_files))

        subset_cram_files = []
        cram_files = list(cram_path.glob("*.cram"))
        console.print(f"[yellow]Found {len(cram_files)} CRAM files in {cram_dir}[/yellow]")

        with ThreadPoolExecutor(max_workers=1) as executor:
            futures = []

            for cram_file in cram_files:
                subset_file = output_path / f"{cram_file.stem}_subset.cram"
                console.print(f"[blue]Subsetting {cram_file.name} to {region}[/blue]")
                futures.append(executor.submit(process_func, str(cram_file), str(subset_file)))

            for future in as_completed(futures):
                subset_cram_files.append(future.result())
                console.print(f"[green]Created subset CRAM: {Path(future.result()).name}[/green]")
                progress.update(task, advance=1)

    console.print(f"[green]Subset CRAM files saved to {output_dir}[/green]")
    return subset_cram_files