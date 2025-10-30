# grid/utils/align_lpa_dir/extract_reads.py
# In[1]: Imports
import pysam
from typing import List, Dict

# In[2]: Function to extract reads from CRAM/BAM
def extract_reads_from_cram(
    cram_file: str,
    reference_fa: str,
    region: str
) -> List[Dict]:
    """
    Extract properly paired reads from CRAM/BAM file using pysam.
    
    Filters for flags: 99, 147, 83, 163 (proper pairs, forward/reverse)
    
    Args:
        cram_file (str): Path to CRAM/BAM file
        reference_fa (str): Reference genome FASTA
        region (str): Genomic region to extract (e.g., 'chr6:160000000-161000000')
        
    Returns:
        List[Dict]: List of read dictionaries with keys:
            - qname: read name
            - pos: alignment position (1-based)
            - cigar: CIGAR string
            - seq: sequence
            - qual: quality string
    """
    reads = []
    
    # Proper pair flags we want
    PROPER_PAIR_FLAGS = {99, 147, 83, 163}
    
    try:
        # Open CRAM/BAM file with reference
        with pysam.AlignmentFile(cram_file, 'rc', reference_filename=reference_fa) as bam:
            # Parse region string (e.g., 'chr6:160000000-161000000')
            if ':' in region:
                chrom, coords = region.split(':')
                start, end = map(int, coords.split('-'))
            else:
                raise ValueError(f"Invalid region format: {region}")
            
            # Fetch reads in region
            for read in bam.fetch(chrom, start, end):
                # Filter for proper pair flags
                if read.flag not in PROPER_PAIR_FLAGS:
                    continue
                
                # Skip unmapped reads
                if read.is_unmapped:
                    continue
                
                # Skip reads without CIGAR
                if read.cigarstring is None:
                    continue
                
                # Get quality string
                qual = pysam.array_to_qualitystring(read.query_qualities)
                
                reads.append({
                    'qname': read.query_name,
                    'pos': read.reference_start + 1,  # Convert to 1-based
                    'cigar': read.cigarstring,
                    'seq': read.query_sequence,
                    'qual': qual
                })
        
        # Sort by qname and position for consistent pairing
        reads.sort(key=lambda x: (x['qname'], x['pos']))
        
        return reads
        
    except Exception as e:
        raise RuntimeError(f"Failed to extract reads from {cram_file}: {e}")