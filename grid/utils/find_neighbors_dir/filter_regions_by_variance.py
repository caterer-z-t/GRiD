# grid/utils/find_neighbors_dir/filter_regions_by_variance.py
# In[1]: Imports
import numpy as np
from rich.console import Console

# In[2]: Console setup
console = Console()

# In[3]: Function to filter regions by variance
def filter_regions_by_variance(data_matrix, regions, sigma2_max=1000.0, frac_r=1.0):
    """
    Filter regions based on variance ratios, similar to C++ code.

    Args:
        data_matrix (np.ndarray): Data matrix of shape [regions x individuals].
        regions (list): List of region identifiers.
        sigma2_max (float): Maximum allowed variance ratio.
        frac_r (float): Fraction of regions to retain based on variance.

    Returns:
        valid_indices (np.ndarray): Indices of regions that pass the filtering.
        variance_ratios (np.ndarray): Calculated variance ratios for all regions.
    """
    console.print(f"[yellow]Filtering regions based on variance...[/yellow]")
    
    R, N = data_matrix.shape
    variance_ratios = []
    
    # Calculate variance ratio for each region
    for r in range(R):
        region_data = data_matrix[r, :]
        # Remove NaN values
        valid_data = region_data[~np.isnan(region_data)]
        
        if len(valid_data) == 0:
            variance_ratios.append(np.inf)
            continue
            
        mean_val = np.mean(valid_data)
        var_val = np.var(valid_data)
        
        # Calculate variance ratio (similar to sigma2ratio in C++)
        if mean_val != 0:
            variance_ratios.append(abs(var_val))  # Use absolute variance as proxy
        else:
            variance_ratios.append(np.inf)
    
    variance_ratios = np.array(variance_ratios)
    
    # Sort and determine thresholds (similar to C++ logic)
    sorted_vars = np.sort(variance_ratios[np.isfinite(variance_ratios)])
    if len(sorted_vars) == 0:
        console.print("[red]Warning: No valid variance ratios found[/red]")
        return np.arange(R), variance_ratios
    
    sigma2_min = sorted_vars[int(len(sorted_vars) * (1 - frac_r))] if frac_r < 1.0 else 0
    
    # Find regions that pass filtering
    valid_indices = []
    extreme_count = 0
    
    for r in range(R):
        var_ratio = variance_ratios[r]
        if np.isfinite(var_ratio) and sigma2_min <= var_ratio <= sigma2_max:
            valid_indices.append(r)
        elif var_ratio > sigma2_max:
            extreme_count += 1
    
    valid_indices = np.array(valid_indices)
    
    console.print(f"[green]✓ Filtering complete[/green]")
    console.print(f"[blue]Kept {len(valid_indices)} out of {R} regions[/blue]")
    console.print(f"[red]Removed {extreme_count} of {R} regions with variance > {sigma2_max}[/red]")
    console.print(f"[yellow]Keeping {len(valid_indices)} regions for analysis[/yellow]")

    return valid_indices, variance_ratios