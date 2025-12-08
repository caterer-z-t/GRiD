# grid/utils/compute_dipcn.py
# In[1]: Imports
from pathlib import Path
import time
from rich.console import Console
from rich.progress import Progress, SpinnerColumn, TextColumn, BarColumn

from .compute_dipcn_dir.load_count_results import load_count_results
from .compute_dipcn_dir.load_neighbor_results import load_neighbor_results
from .compute_dipcn_dir.validate_sample_overlap import validate_sample_overlap
from .compute_dipcn_dir.compute_diploid_cn import compute_diploid_cn_for_exon
from .compute_dipcn_dir.write_dipcn_output import write_dipcn_output

# Exon types to process
EXON_TYPES = ['1B_KIV3', '1B_notKIV3', '1B', '1A']

# In[2]: Main function to compute diploid copy numbers
def compute_dipcn_pipeline(
    count_file: Path,
    neighbor_file: Path,
    output_prefix: Path,
    n_neighbors: int,
    console: Console = None
):
    """
    Main function to compute diploid copy numbers for all exon types.
    
    Args:
        count_file (Path): Path to realignment count file
        neighbor_file (Path): Path to neighbor results file (can be gzipped)
        output_prefix (Path): Output file prefix
        n_neighbors (int): Number of top neighbors to use
        console (Console): Rich console for output (optional)
        
    Returns:
        None
    """
    if console is None:
        console = Console()

    console.print(f"[blue]Starting diploid copy number computation...[/blue]")

    output_prefix = Path(output_prefix).expanduser()
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    
    start_time = time.time()
    
    # Step 1: Load count results
    console.rule("[bold blue]Step 1: Load Count Results")
    counts = load_count_results(count_file)
    console.print(f"[green]✓ Loaded counts for {len(counts)} samples[/green]")
    
    # Step 2: Load neighbor results
    console.rule("[bold blue]Step 2: Load Neighbor Results")
    neighbors = load_neighbor_results(neighbor_file)
    console.print(f"[green]✓ Loaded neighbor info for {len(neighbors)} samples[/green]")
    
    # Step 3: Validate sample overlap
    console.rule("[bold blue]Step 3: Validate Sample Overlap")
    overlap_count, overlap_samples = validate_sample_overlap(counts, neighbors, console)
    
    if overlap_count == 0:
        console.print("[red]✗ No overlapping samples found! Cannot proceed.[/red]")
        raise ValueError("No overlapping sample IDs between count and neighbor files")
    
    console.print(f"[green]✓ Found {overlap_count} overlapping samples[/green]")
    
    # Step 4: Create output directory
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    
    # Step 5: Process each exon type
    console.rule("[bold blue]Step 4: Compute Diploid Copy Numbers")
    
    results_summary = {}
    
    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        console=console,
    ) as progress:
        task = progress.add_task(
            f"[cyan]Processing exon types...",
            total=len(EXON_TYPES)
        )
        
        for exon_type in EXON_TYPES:
            # Compute diploid copy numbers for this exon type
            results = compute_diploid_cn_for_exon(
                counts=counts,
                neighbors=neighbors,
                exon_type=exon_type,
                n_neighbors=n_neighbors
            )
            
            # Write output
            output_file = f"{output_prefix}.exon{exon_type}.dipCN.txt"
            write_dipcn_output(results, output_file)
            
            results_summary[exon_type] = {
                'count': len(results),
                'file': output_file
            }
            
            progress.advance(task)
    
    # Step 6: Summary
    console.rule("[bold blue]Summary")
    elapsed = time.time() - start_time
    
    console.print(f"[bold green]✓ Processed {len(EXON_TYPES)} exon types[/bold green]")
    console.print("")
    
    for exon_type, info in results_summary.items():
        console.print(f"[blue]{exon_type}:[/blue] {info['count']} samples → {info['file']}")
    
    console.print("")
    console.print(f"[blue]Time elapsed: {elapsed:.2f} seconds[/blue]")
    console.rule("[bold green]✓ Diploid copy number computation complete")