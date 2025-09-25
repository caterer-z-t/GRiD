# In[0]: Imports
import argparse
import pysam
from pathlib import Path

# In[1]: Functions
def format_region(args) -> str:
    """
    Format region string from chromosome, start, and end.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command-line arguments.
    Returns
    -------
    str
        Formatted region string ("chr:start-end").
    """

    if args.chromosome and not args.chromosome.startswith("chr"):
        args.chromosome = f"chr{args.chromosome}"

    return f"{args.chromosome}:{args.start}-{args.end}"

def arg_check(args):

    if not Path(args.input).exists():
        raise FileNotFoundError(f"Input file {args.input} does not exist.")
    
    if args.input.endswith(".cram") and args.reference is None:
        raise ValueError("Reference FASTA must be provided for CRAM files.")
    
    if not Path(args.input + ".crai").exists() and args.input.endswith(".cram"):
        raise FileNotFoundError(f"CRAM index file {args.input}.crai does not exist.")
    
    if not Path(args.input + ".bai").exists() and args.input.endswith(".bam"):
        raise FileNotFoundError(f"BAM index file {args.input}.bai does not exist.")
    
    if args.region is None and (args.start is None or args.end is None or args.chromosome is None):
        raise ValueError("Either --region or all of --start, --end, and --chromosome must be provided.")


def argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run mosdepth on one CRAM file for a single region."
    )
    parser.add_argument(
        "-i", "--input", required=True, help="Input CRAM file (must be indexed)"
    )
    parser.add_argument(
        "-o", "--output", required=True, help="Output TSV file"
    )
    parser.add_argument(
        "-R", "--reference", required=True, help="Reference FASTA (required for CRAM)"
    )
    parser.add_argument(
        "-r", "--region", required=False, help="Genomic region (e.g., chr6:135000000-135100000)"
    )
    parser.add_argument(
        "-c", "--chromosome", required=False, help="Chromosome name (e.g., chr6)"
    )
    parser.add_argument(
        "-s", "--start", required=False, type=int, help="Start position (1-based)"
    )
    parser.add_argument(
        "-e", "--end", required=False, type=int, help="End position (1-based)"
    )
    return parser

def ensure_index(bam_or_cram: str, reference: str | None = None) -> str:
    """
    Ensure that a BAM/CRAM file has an index. If missing, create it.

    Parameters
    ----------
    bam_or_cram : str
        Path to BAM or CRAM file.
    reference : str | None
        Path to reference FASTA (only needed for CRAM indexing if required).

    Returns
    -------
    str
        Path to the index file (.bai for BAM, .crai for CRAM).
    """
    path = Path(bam_or_cram)

    if not path.exists():
        raise FileNotFoundError(f"Input file {bam_or_cram} does not exist.")

    if bam_or_cram.endswith(".bam"):
        index_path = path.with_suffix(path.suffix + ".bai")
    elif bam_or_cram.endswith(".cram"):
        index_path = path.with_suffix(path.suffix + ".crai")
    else:
        raise ValueError("Input file must be BAM or CRAM.")

    if not index_path.exists():
        pysam.index(str(path), reference=reference)
    
    return str(index_path)
