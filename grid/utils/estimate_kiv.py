# grid/utils/estimate_kiv.py
# In[1]: Imports
import time
import pandas as pd
from rich.console import Console
from rich.table import Table
from pathlib import Path

from .estimate_kiv_dir.load_dipcn_file import load_dipcn_file
from .estimate_kiv_dir.merge_exon_data import merge_exon_data
from .estimate_kiv_dir.compute_estimates import compute_kiv2_estimates
from .estimate_kiv_dir.compute_summary_stats import compute_summary_stats

# In[2]: Main Function
def compute_kiv_estimates(
    output: str,
    exon1a_file: str,
    exon1b_file: str,
    console: Console = None
) -> pd.DataFrame:
    """
    Compute KIV2 copy number estimates from exon1A and exon1B diploid copy numbers.
    
    Formula:
        diploid_estimate = 34.9 × exon1A + 5.2 × exon1B - 1
        haploid_estimate = diploid_estimate / 2
    
    Args:
        exon1a_file (str): Path to exon1A diploid CN file
        exon1b_file (str): Path to exon1B diploid CN file
        console (Console): Rich console for output (optional)
        
    Returns:
        pd.DataFrame: DataFrame with columns [exon1A, exon1B, dip_estimate, estimate]
        
    Raises:
        ValueError: If no overlapping samples between files
    """

        
    # Handle path operations in CLI
    output_path = Path(output).expanduser().resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    if console is None:
        console = Console()
    
    console.print(f"[blue]Starting KIV2 copy number estimation...[/blue]")

    start_time = time.time()
    
    # Step 1: Load exon1A data
    console.rule("[bold blue]Step 1: Load Exon1A Diploid Copy Numbers")
    exon1a = load_dipcn_file(exon1a_file, 'exon1A')
    console.print(f"[green]✓ Loaded {len(exon1a)} samples[/green]")
    
    # Step 2: Load exon1B data
    console.rule("[bold blue]Step 2: Load Exon1B Diploid Copy Numbers")
    exon1b = load_dipcn_file(exon1b_file, 'exon1B')
    console.print(f"[green]✓ Loaded {len(exon1b)} samples[/green]")
    
    # Step 3: Merge data
    console.rule("[bold blue]Step 3: Merge Exon Data")
    combined = merge_exon_data(exon1a, exon1b)
    
    if len(combined) == 0:
        console.print("[red]✗ No overlapping samples found![/red]")
        raise ValueError("No overlapping samples between exon1A and exon1B files")
    
    console.print(f"[green]✓ Found {len(combined)} samples with both exon1A and exon1B[/green]")
    
    # Step 4: Compute estimates
    console.rule("[bold blue]Step 4: Compute KIV2 Copy Number Estimates")
    results = compute_kiv2_estimates(combined)
    console.print(f"[green]✓ Computed estimates for {len(results)} samples[/green]")
    
    # Step 5: Display summary statistics
    console.rule("[bold blue]Summary Statistics")
    stats = compute_summary_stats(results)
    
    # Create a nice table
    table = Table(show_header=True, header_style="bold magenta")
    table.add_column("Metric", style="cyan")
    table.add_column("Diploid Estimate", justify="right")
    table.add_column("Haploid Estimate", justify="right")
    
    table.add_row("Mean", f"{stats['dip_mean']:.2f}", f"{stats['hap_mean']:.2f}")
    table.add_row("Median", f"{stats['dip_median']:.2f}", f"{stats['hap_median']:.2f}")
    table.add_row("Std Dev", f"{stats['dip_std']:.2f}", f"{stats['hap_std']:.2f}")
    table.add_row("Min", f"{stats['dip_min']:.2f}", f"{stats['hap_min']:.2f}")
    table.add_row("Max", f"{stats['dip_max']:.2f}", f"{stats['hap_max']:.2f}")
    
    console.print(table)
    
    elapsed = time.time() - start_time
    console.print(f"\n[blue]Time elapsed: {elapsed:.2f} seconds[/blue]")
    console.rule("[bold green]✓ KIV2 estimation complete")


    # Save results (CLI handles file writing)
    console.rule("[bold blue]Step 5: Write Output")
    console.print(f"[yellow]Writing results to {output_path}[/yellow]")
    
    if format == 'csv':
        results.to_csv(output_path, index=True)
    else:  # tsv or txt
        results.to_csv(output_path, index=True, sep='\t')
    
    console.print(f"[green]✓ Wrote {len(results)} samples to {output_path}[/green]")
    
    # Display first few rows
    console.print("\n[bold]First 10 samples:[/bold]")
    
    table = Table(show_header=True, header_style="bold cyan")
    table.add_column("Sample ID", style="cyan")
    table.add_column("exon1A", justify="right")
    table.add_column("exon1B", justify="right")
    table.add_column("dip_estimate", justify="right")
    table.add_column("estimate", justify="right")
    
    for idx, row in results.head(10).iterrows():
        table.add_row(
            str(idx),
            f"{row['exon1A']:.4f}",
            f"{row['exon1B']:.4f}",
            f"{row['dip_estimate']:.2f}",
            f"{row['estimate']:.2f}"
        )
    
    console.print(table)
    console.print(f"\n[green]✓ KIV2 estimation completed successfully![/green]")
    return results