# grid/utils/estimate_kiv_dir/compute_summary_stats.py
# In[1]: Imports
import pandas as pd
from typing import Dict

# In[2]: Compute Summary Statistics
def compute_summary_stats(results: pd.DataFrame) -> Dict[str, float]:
    """
    Compute summary statistics for diploid and haploid estimates.
    
    Args:
        results (pd.DataFrame): DataFrame with columns ['dip_estimate', 'estimate']
        
    Returns:
        Dict[str, float]: Dictionary with summary statistics
    """
    stats = {
        # Diploid statistics
        'dip_mean': results['dip_estimate'].mean(),
        'dip_median': results['dip_estimate'].median(),
        'dip_std': results['dip_estimate'].std(),
        'dip_min': results['dip_estimate'].min(),
        'dip_max': results['dip_estimate'].max(),
        
        # Haploid statistics
        'hap_mean': results['estimate'].mean(),
        'hap_median': results['estimate'].median(),
        'hap_std': results['estimate'].std(),
        'hap_min': results['estimate'].min(),
        'hap_max': results['estimate'].max(),
    }
    
    return stats