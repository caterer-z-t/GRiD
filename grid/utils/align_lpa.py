# grid/utils/align_lpa.py
# In[1]: Imports
from pathlib import Path
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from threading import Lock
from functools import partial

from rich.console import Console
from rich.progress import Progress, SpinnerColumn, TextColumn, BarColumn, TaskProgressColumn

from .align_lpa_dir.find_cram_files import find_cram_files
from .align_lpa_dir.load_positions import load_positions
from .align_lpa_dir.load_reference import load_reference
from .align_lpa_dir.process_cram_file import process_cram_file
from .align_lpa_dir.validate_inputs import validate_inputs

# In[2]: Set up console for rich output
console = Console()

# In[3]: Main function to run LPA realignments
def run_lpa_realignments(
    cram_dir,
    reference_fa,
    lpa_ref_fasta,
    positions_file,
    genome_build,
    chrom,
    start,
    end,
    output_file,
    threads: int = 1
):
    """
    Main function to run LPA realignment for all CRAM files in a directory.
    
    UPDATED: Now supports parallel processing with multiple threads.
    
    Args:
        cram_dir (str): Directory containing CRAM/BAM files
        reference_fa (str): Reference genome FASTA for CRAM files
        lpa_ref_fasta (str): LPA reference FASTA for realignment
        positions_file (str): File with hardcoded repeat positions
        genome_build (str): Genome build ('hg19', 'hg37', or 'hg38')
        chrom (str): Chromosome (e.g., 'chr6' or '6')
        start (int): Start position (0-based)
        end (int): End position (0-based)
        output_file (str): Output TSV file path
        threads (int): Number of threads for parallel processing (default: 1)
        
    Returns:
        None
    """
    start_time = time.time()
    
    # Step 1: Validate all inputs
    console.rule("[bold blue]Step 1: Validate Inputs")
    validate_inputs(
        cram_dir=cram_dir,
        reference_fa=reference_fa,
        lpa_ref_fasta=lpa_ref_fasta,
        positions_file=positions_file
    )
    console.print("[green]✓ All input files validated[/green]")
    
    # Step 2: Find CRAM/BAM files
    console.rule("[bold blue]Step 2: Find CRAM/BAM Files")
    cram_files = find_cram_files(cram_dir)
    console.print(f"[green]✓ Found {len(cram_files)} CRAM/BAM files[/green]")
    
    if len(cram_files) == 0:
        console.print("[red]✗ No CRAM/BAM files found[/red]")
        raise FileNotFoundError(f"No CRAM/BAM files in {cram_dir}")
    
    # Step 3: Load positions
    console.rule("[bold blue]Step 3: Load Repeat Positions")
    # Treat hg37 as hg19
    build = 'hg19' if genome_build == 'hg37' else genome_build
    starts, ref_offset = load_positions(positions_file, build)
    console.print(f"[green]✓ Loaded {len(starts)} repeat positions[/green]")
    console.print(f"[blue]  Reference offset: {ref_offset}[/blue]")
    console.print(f"[blue]  Repeat starts: {starts[:3]}...[/blue]")
    
    # Step 4: Load reference sequence
    console.rule("[bold blue]Step 4: Load LPA Reference Sequence")
    ref_seq = load_reference(lpa_ref_fasta)
    console.print(f"[green]✓ Loaded reference: {len(ref_seq):,} bp[/green]")
    
    # Step 5: Create output directory and file
    output_path = Path(output_file).expanduser().resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Remove existing output file if present
    if output_path.exists():
        output_path.unlink()
        console.print(f"[yellow]Removed existing output file[/yellow]")
    
    # Step 6: Process all CRAM files (with multithreading)
    console.rule("[bold blue]Step 5: Process CRAM Files")
    
    # Determine region string (with or without 'chr' prefix)
    region = f"{chrom}:{start}-{end}"
    
    # Display threading information
    console.print(f"[yellow]Using {threads} thread(s) for parallel processing[/yellow]")
    
    # Create lock for thread-safe file writing
    write_lock = Lock()
    
    # Create partial function with fixed parameters
    process_func = partial(
        process_cram_file,
        reference_fa=reference_fa,
        region=region,
        ref_seq=ref_seq,
        starts=starts,
        ref_offset=ref_offset,
        output_file=str(output_path)
    )
    
    # Process CRAMs with threading
    processed = 0
    failed = 0
    total_reads = 0
    failed_samples = []
    
    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        TaskProgressColumn(),
        console=console,
    ) as progress:
        task = progress.add_task(
            f"[cyan]Processing files...",
            total=len(cram_files)
        )
        
        with ThreadPoolExecutor(max_workers=threads) as executor:
            # Submit all jobs
            future_to_cram = {
                executor.submit(process_func, cram_file): cram_file
                for cram_file in cram_files
            }
            
            # Collect results as they complete
            for future in as_completed(future_to_cram):
                cram_file = future_to_cram[future]
                
                try:
                    result = future.result()
                    
                    if result['success']:
                        processed += 1
                        total_reads += result.get('read_count', 0)
                    else:
                        failed += 1
                        failed_samples.append(cram_file.name)
                        error_msg = result.get('error', 'Unknown error')
                        progress.console.print(
                            f"[red]✗ {cram_file.name}: {error_msg}[/red]"
                        )
                        
                except Exception as e:
                    failed += 1
                    failed_samples.append(cram_file.name)
                    progress.console.print(
                        f"[red]✗ Failed to process {cram_file.name}: {e}[/red]"
                    )
                
                progress.update(task, advance=1)
    
    # Step 7: Summary
    console.rule("[bold blue]Summary")
    elapsed = time.time() - start_time
    
    console.print(f"[green]✓ Successfully processed: {processed} files[/green]")
    if failed > 0:
        console.print(f"[yellow]⚠ Failed/skipped: {failed} files[/yellow]")
        if len(failed_samples) <= 10:
            for sample in failed_samples:
                console.print(f"    - {sample}")
        else:
            for sample in failed_samples[:5]:
                console.print(f"    - {sample}")
            console.print(f"    ... and {len(failed_samples) - 5} more")
    
    console.print(f"[blue]Total filtered reads: {total_reads:,}[/blue]")
    console.print(f"[blue]Time elapsed: {elapsed:.2f} seconds[/blue]")
    
    if threads > 1:
        avg_time = elapsed / len(cram_files)
        console.print(f"[blue]Average time per file: {avg_time:.2f} seconds[/blue]")
    
    console.print(f"[bold green]✓ Output written to: {output_path}[/bold green]")
    
    # Show sample of output
    if output_path.exists():
        with open(output_path, 'r') as f:
            lines = [next(f, None) for _ in range(5)]
            lines = [line for line in lines if line]  # Remove None values
        
        if lines:
            console.print("\n[bold]First few lines of output:[/bold]")
            for line in lines:
                console.print(f"  {line.strip()}")
