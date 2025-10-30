# grid/utils/compute_dipcn_dir/normalize_sample_id.py
# In[1]: Function to normalize sample IDs
def normalize_sample_id(sample_id: str) -> str:
    """
    Normalize sample IDs by removing common file extensions and suffixes.
    
    Examples:
        NWD278973.b38.irc.v1_subset → NWD278973
        NWD278973.cram → NWD278973
        NWD278973 → NWD278973
        
    Args:
        sample_id (str): Original sample ID
        
    Returns:
        str: Normalized sample ID
    """
    # Strip whitespace first
    sample_id = sample_id.strip()
    
    # Remove .b38.irc.v1_subset pattern (before other extensions)
    if '.b38.irc.v1_subset' in sample_id:
        sample_id = sample_id.replace('.b38.irc.v1_subset', '')
    
    # Remove common CRAM/BAM file patterns
    if sample_id.endswith('.cram'):
        sample_id = sample_id[:-5]
    elif sample_id.endswith('.bam'):
        sample_id = sample_id[:-4]
    
    return sample_id.strip()