# grid/utils/run_mosdepth.py
# In[1]: Imports
import sys
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from threading import Lock
from functools import partial
import csv

from rich.console import Console
from rich.progress import Progress, SpinnerColumn, TextColumn, BarColumn, TaskProgressColumn

from .mosdepth_dir.check_mosdepth_avail import check_mosdepth_available
from .helper_dir.find_all_cram_files import find_cram_files
from .mosdepth_dir.run_mosdepth_single_cram import run_mosdepth_single_cram
from .mosdepth_dir.write_coverage_result import write_coverage_result
from .helper_dir.display_results import print_batch_summary

# In[2]: Setup Rich Console
console = Console()

# In[3]: Main Function to Run Mosdepth
def run_mosdepth(
    cram_dir: str,
    output_file: str,
    ref_fasta: str,
    chrom: str,
    start: int,
    end: int,
    region_name: str,
    work_dir: str,
    by: int,
    fast_mode: bool,
    threads: int = 1
):
    """
    Run mosdepth on all CRAMs in a directory using parallel processing.
    
    Args:
        cram_dir: Directory containing CRAM files
        output_file: Path to output TSV file
        ref_fasta: Path to reference genome FASTA
        chrom: Chromosome name
        start: Start position
        end: End position
        work_dir: Working directory for intermediate files
        by: Bin size for mosdepth
        fast_mode: Whether to use fast mode
        threads: Number of threads for parallel processing
    """
    # Check mosdepth availability
    try:
        check_mosdepth_available()
    except RuntimeError as e:
        console.print(f"[red]✗ {str(e)}[/red]")
        sys.exit(1)
    
    # Setup directories
    work_path = Path(work_dir).expanduser()
    work_path.mkdir(parents=True, exist_ok=True)
    output_path = Path(output_file).expanduser()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Find CRAM files
    cram_files = find_cram_files(cram_dir)
    if not cram_files:
        console.print(f"[red]✗ No CRAM files found in {cram_dir}[/red]")
        sys.exit(1)
    
    # Setup output file with header
    with open(output_path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["Sample", f"{chrom}:{start}-{end}"])
    
    # Display initial information
    console.print(f"[yellow]Detected {len(cram_files)} CRAM files to process[/yellow]")
    console.print(f"[yellow]Using {threads} thread(s) for parallel processing[/yellow]")
    
    # Create lock for thread-safe file writing
    write_lock = Lock()
    
    # Create partial function with fixed parameters
    process_func = partial(
        run_mosdepth_single_cram,
        ref_fasta=ref_fasta,
        work_dir=work_path,
        chrom=chrom,
        start=start,
        end=end,
        region_name=region_name,
        by=by,
        fast_mode=fast_mode
    )
    
    # Process CRAMs with threading and progress tracking
    failed = []
    
    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        TaskProgressColumn(),
        console=console
    ) as progress:
        
        task = progress.add_task(
            "[bold green]Running mosdepth", 
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
                sample_name, coverage = future.result()
                
                if coverage != "Error":
                    write_coverage_result(
                        output_path, 
                        sample_name, 
                        coverage, 
                        write_lock,
                        progress.console
                    )
                else:
                    progress.console.print(f"[red]✗ {sample_name} failed[/red]")
                    failed.append(sample_name)
                
                progress.update(task, advance=1)
    
    # Print summary
    print_batch_summary(len(cram_files), failed, operation="processed")
    console.print(f"[bold blue]Results written to {output_path}[/bold blue]")