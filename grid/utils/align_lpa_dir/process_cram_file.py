# grid/utils/align_lpa_dir/process_cram_file.py
# In[1]: Imports
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, List

from .extract_reads import extract_reads_from_cram
from .classify_reads import classify_read_pairs

# In[2]: Function to process CRAM/BAM file
def process_cram_file(
    cram_file: Path,
    reference_fa: str,
    region: str,
    ref_seq: str,
    starts: List[int],
    ref_offset: int,
    output_file: str
) -> Dict:
    """
    Process a single CRAM/BAM file to classify LPA variants.
    
    Args:
        cram_file (Path): Path to CRAM/BAM file
        reference_fa (str): Reference genome FASTA
        region (str): Genomic region (e.g., 'chr6:160000000-161000000')
        ref_seq (str): LPA reference sequence
        starts (List[int]): List of 7 repeat start positions
        ref_offset (int): Reference offset position
        output_file (str): Output file to append results
        
    Returns:
        Dict: Result dictionary with 'success', 'sample_id', 'read_count', 'counts'
    """
    sample_id = cram_file.stem  # Remove .cram or .bam extension
    
    try:
        # Step 1: Extract reads from CRAM using samtools
        reads = extract_reads_from_cram(
            cram_file=str(cram_file),
            reference_fa=reference_fa,
            region=region
        )
        
        if len(reads) == 0:
            # No reads in region - write zeros
            with open(output_file, 'a') as out:
                out.write(f"{sample_id}\t0\t0\t0\t0\n")
            
            return {
                'success': True,
                'sample_id': sample_id,
                'read_count': 0,
                'counts': {'1B_KIV3': 0, '1B_KIV2': 0, '1B_tied': 0, '1A': 0}
            }
        
        # Step 2: Classify read pairs
        counts = classify_read_pairs(
            reads=reads,
            ref_seq=ref_seq,
            starts=starts,
            ref_offset=ref_offset
        )
        
        # Step 3: Write results to output file
        with open(output_file, 'a') as out:
            out.write(
                f"{sample_id}\t{counts['1B_KIV3']}\t{counts['1B_KIV2']}\t"
                f"{counts['1B_tied']}\t{counts['1A']}\n"
            )
        
        return {
            'success': True,
            'sample_id': sample_id,
            'read_count': len(reads),
            'counts': counts
        }
        
    except Exception as e:
        # Write zeros on error
        with open(output_file, 'a') as out:
            out.write(f"{sample_id}\t0\t0\t0\t0\n")
        
        return {
            'success': False,
            'sample_id': sample_id,
            'error': str(e)
        }