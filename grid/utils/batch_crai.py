# GRiD/grid/utils/batch_crai.py
# In[1]: Imports
from pathlib import Path
import pysam
from concurrent.futures import ThreadPoolExecutor, as_completed
from threading import Lock
from functools import partial
from rich.console import Console
from rich.progress import Progress, SpinnerColumn, TextColumn, BarColumn, TaskProgressColumn

from .ensure_crai import ensure_crai

# In[2]: Function to Batch Create CRAI Indexes
def batch_crai(cram_dir, reference, console, threads=1):
    """
    Batch create CRAI indexes for all CRAM files in a directory.

    Args:
        cram_dir: Directory containing CRAM files.
        reference: Reference genome FASTA for indexing.
        console: Console object for logging.
        threads: Number of threads for parallel processing. Default is 1.

    Returns:
        Path to CRAI files.
    """
    cram_path = Path(cram_dir)

    if not cram_path.is_dir():
        raise NotADirectoryError(f"Provided path is not a directory: {cram_dir}")

    process_func = partial(
        ensure_crai,
        reference=reference
        )

    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        TaskProgressColumn(),
        console=console
    ) as progress:

        cram_list = list(cram_path.glob("*.cram"))
        task = progress.add_task("[bold green] Creating CRAI indexes...", total=len(cram_list))

        with ThreadPoolExecutor(max_workers=threads) as executor:
            futures = []

            for cram_file in cram_list:
                crai_file = cram_file.with_suffix(cram_file.suffix + ".crai")

                if crai_file.exists():
                    progress.update(task, advance=1)
                    continue

                futures.append(executor.submit(process_func, str(cram_file)))

            for future in as_completed(futures):
                cram_file = future.result()
                crai_file = cram_file.with_suffix(cram_file.suffix + ".crai")

                progress.console.log(f"[green]Created CRAI: {crai_file.name}[/green]")
                progress.update(task, advance=1)

    return Path(cram_dir)
