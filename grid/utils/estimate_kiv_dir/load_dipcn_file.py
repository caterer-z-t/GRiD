# grid/utils/estimate_kiv_dir/load_dipcn_file.py
# In[1]: Imports
import pandas as pd

# In[2]: Load Diploid Copy Number File
def load_dipcn_file(file_path: str, exon_name: str) -> pd.DataFrame:
    """
    Load diploid copy number file.
    
    Expected format:
        ID    dipCN
        sample1    15.234
        sample2    18.456
    
    Args:
        file_path (str): Path to diploid CN file
        exon_name (str): Name to use for the column (e.g., 'exon1A', 'exon1B')
        
    Returns:
        pd.DataFrame: DataFrame with index=ID and column=exon_name
        
    Raises:
        FileNotFoundError: If file doesn't exist
        ValueError: If file format is invalid
    """
    try:
        # Try tab-separated first (our output format)
        df = pd.read_csv(file_path, sep='\t', index_col=0)
    except Exception as e:
        try:
            # Try space-separated (alternative format)
            df = pd.read_csv(file_path, sep=' ', index_col=0)
        except Exception as e2:
            raise ValueError(f"Failed to load {file_path}: {e}") from e2
    
    # Rename column to exon name
    df.rename(columns={'dipCN': exon_name}, inplace=True)
    
    return df