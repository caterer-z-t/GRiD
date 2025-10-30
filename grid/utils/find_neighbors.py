# grid/utils/find_neighbors.py
# In[1]: Imports
from pathlib import Path
import sys
import time
from rich.console import Console
from rich.progress import Progress, SpinnerColumn, TextColumn, BarColumn, TaskProgressColumn

from .find_neighbors_dir.read_normalized_data import read_normalized_data
from .find_neighbors_dir.filter_regions_by_variance import filter_regions_by_variance
from .find_neighbors_dir.find_neighbors_sklearn import find_neighbors_sklearn
from .find_neighbors_dir.save_neighbors import save_neighbors

# In[2]: Console setup
console = Console()

# In[3]: Main function to find neighbors
def find_neighbors(input_file, output_file, zmax, n_neighbors, sigma2_max):
    """
    Main function to find nearest neighbors in normalized depth data.

    Args:
        input_file (str): Path to gzipped normalized depth file.
        output_file (str): Path to output file (will be converted to .zMax{zmax}.txt.gz).
        zmax (float): Maximum z-score clipping threshold.
        n_neighbors (int): Number of nearest neighbors to find.
        sigma2_max (float): Maximum variance ratio for region filtering.

    Returns:
        None
    """
    start_time = time.time()
    output_path = Path(output_file).expanduser().resolve()
    
    # Extract output directory and prefix
    output_dir = output_path.parent
    output_prefix = str(output_dir / output_path.stem)

    console.rule("[bold blue]Step 1: Read Normalized Data")
    individuals, regions, data_matrix, scales = read_normalized_data(input_file)
    console.print(f"[green]✓ Loaded {len(individuals)} individuals × {len(regions)} regions[/green]")

    console.rule("[bold blue]Step 2: Filter Regions by Variance")
    valid_indices, variance_ratios = filter_regions_by_variance(
        data_matrix, regions, sigma2_max=sigma2_max
    )
    R_use = len(valid_indices)
    console.print(f"[green]✓ Retained {R_use:,} regions below variance ratio {sigma2_max}[/green]")

    console.rule("[bold blue]Step 3: Compute Nearest Neighbors")
    
    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        TaskProgressColumn(),
        console=console,
    ) as progress:
        task = progress.add_task("[cyan]Finding neighbors...", total=1)
        
        # Use filtered data matrix (regions x individuals)
        filtered_data = data_matrix[valid_indices, :]
        
        # Find neighbors for all individuals at once
        neighbors = find_neighbors_sklearn(
            filtered_data,
            individuals,
            n_neighbors=n_neighbors,
            zmax=zmax
        )
        
        progress.advance(task)

    console.print(f"[green]✓ Neighbor search complete for all {len(individuals)} samples[/green]")

    console.rule("[bold blue]Step 4: Save Output")
    save_neighbors(neighbors, scales, output_prefix, zmax, R_use)
    console.print(f"[bold green]✓ Output written to {output_prefix}.zMax{zmax:.1f}.txt.gz[/bold green]")

    elapsed = time.time() - start_time
    console.rule(f"[bold blue]✓ Done in {elapsed:.2f} seconds[/bold blue]")