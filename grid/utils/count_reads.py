# grid/utils/count_reads_for_directory.py
# In[1]: Imports
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from threading import Lock
from functools import partial

from .helper_dir.find_all_cram_files import find_cram_files
from .helper_dir.setup_output_file import setup_output_file
from .helper_dir.write_result_to_file import write_result_to_file
from .helper_dir.load_flags_from_yaml import load_flags
from .count_reads_dir.process_single_cram import process_single_cram

from rich.console import Console
from rich.progress import Progress, SpinnerColumn, TextColumn, BarColumn, TaskProgressColumn

# In[2]: Setup Rich Console
console = Console()

# In[3]: Main Function to Count Reads in Directory
def count_reads(
    cram_dir: str,
    output_file: str,
    ref_fasta: str,
    chrom: str,
    start: int,
    end: int,
    config_file: str,
    threads: int = 1
):
    """
    Count reads for all CRAMs in a directory using parallel processing.
    
    Args:
        cram_dir: Directory containing CRAM files
        output_file: Path to output TSV file
        ref_fasta: Path to reference genome FASTA
        chrom: Chromosome name
        start: Start position
        end: End position
        config_file: Path to YAML config file
        threads: Number of threads for parallel processing
    """
    # Load configuration
    proper_flags = load_flags(config_file, "read-count")
    
    # Find all CRAM files
    cram_files = find_cram_files(cram_dir)
    
    # Setup output file with header
    output_path = setup_output_file(output_file, chrom, start, end)
    
    # Display initial information
    console.print(f"[yellow]Detected {len(cram_files)} CRAM files in {cram_dir}[/yellow]")
    console.print(f"[yellow]Using {threads} thread(s) for parallel processing[/yellow]")
    
    # Create lock for thread-safe file writing
    write_lock = Lock()
    
    # Create a partial function with fixed parameters
    process_func = partial(
        process_single_cram,
        ref_fasta=ref_fasta,
        chrom=chrom,
        start=start,
        end=end,
        proper_flags=proper_flags
    )
    
    # Process CRAMs with threading and progress tracking
    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        TaskProgressColumn(),
        console=console
    ) as progress:
        
        task = progress.add_task(
            "[bold green]Processing CRAM files", 
            total=len(cram_files)
        )
        
        with ThreadPoolExecutor(max_workers=threads) as executor:
            # Submit all jobs
            future_to_cram = {
                executor.submit(process_func, cram): cram 
                for cram in cram_files
            }
            
            # Collect and write results as they complete
            for future in as_completed(future_to_cram):
                basename, count = future.result()
                
                # Write result to file (thread-safe)
                write_result_to_file(output_path, basename, count, write_lock, progress.console)
                
                # Update progress bar
                progress.update(task, advance=1)
    
    console.print(
        f"[bold blue]VNTR read counting completed. "
        f"Results written to {output_path}[/bold blue]"
    )