# grid/utils/estimate_kiv_dir/merge_exon_data.py
# In[1]: Imports
import pandas as pd

# In[2]: Merge Exon Data
def merge_exon_data(exon1a: pd.DataFrame, exon1b: pd.DataFrame) -> pd.DataFrame:
    """
    Merge exon1A and exon1B data on sample IDs (inner join).
    
    Args:
        exon1a (pd.DataFrame): DataFrame with exon1A data
        exon1b (pd.DataFrame): DataFrame with exon1B data
        
    Returns:
        pd.DataFrame: Combined DataFrame with both exon1A and exon1B columns
    """
    # Inner join - only keep samples present in both files
    combined = pd.concat([exon1a, exon1b], join='inner', axis=1)
    
    return combined