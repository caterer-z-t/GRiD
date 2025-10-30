# grid/utils/estimate_kiv_dir/compute_estimates.py
# In[1]: Imports
import pandas as pd

# Empirically determined coefficients
COEF_EXON1A = 34.9
COEF_EXON1B = 5.2
INTERCEPT = -1

# In[2]: Compute KIV2 Estimates
def compute_kiv2_estimates(combined: pd.DataFrame) -> pd.DataFrame:
    """
    Compute KIV2 copy number estimates using the empirical formula.
    
    Formula:
        diploid_estimate = 34.9 × exon1A + 5.2 × exon1B - 1
        haploid_estimate = diploid_estimate / 2
    
    Args:
        combined (pd.DataFrame): DataFrame with columns ['exon1A', 'exon1B']
        
    Returns:
        pd.DataFrame: DataFrame with added columns ['dip_estimate', 'estimate']
    """
    # Create a copy to avoid modifying original
    results = combined.copy()
    
    # Compute diploid estimate
    results['dip_estimate'] = (
        COEF_EXON1A * results['exon1A'] + 
        COEF_EXON1B * results['exon1B'] + 
        INTERCEPT
    )
    
    # Compute haploid estimate
    results['estimate'] = results['dip_estimate'] / 2
    
    return results