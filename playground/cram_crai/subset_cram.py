# In[1]: Imports
import sys
import pysam
import os
import multiprocessing
from typing import Optional, Dict
import argparse

# In[2]: Subset CRAM/BAM Function

# In[2.A]: Helper: Determine output mode
def determine_mode(output_file: str) -> str:
    ext = os.path.splitext(output_file)[1].lower()
    if ext == ".bam":
        return "wb"
    elif ext == ".cram":
        return "wc"
    else:
        raise ValueError("Output file must have .bam or .cram extension")


# In[2.B]: Helper: Build open kwargs for pysam
def build_open_kwargs(reference_fasta: Optional[str] = None,
                      input_crai: Optional[str] = None) -> dict:
    kwargs = {}
    if reference_fasta:
        kwargs["reference_filename"] = reference_fasta
    if input_crai:
        kwargs["index_filename"] = input_crai
    return kwargs


# In[2.C]: Helper: Open input/output files
def open_files(input_file: str, output_file: str, out_mode: str,
               open_kwargs: dict, threads: int = 1):
    infile = pysam.AlignmentFile(input_file, "rc", **open_kwargs, threads=threads)
    outfile = pysam.AlignmentFile(output_file, out_mode, header=infile.header, threads=threads)
    return infile, outfile


# In[3]: Core: Subset reads
def write_subset(infile, outfile, region: Dict[str, Optional[int]]):
    if region["start"] is not None and region["end"] is not None:
        for read in infile.fetch(region["chrom"], region["start"], region["end"]):
            outfile.write(read)
    else:
        # Fetch full chromosome if no start/end provided
        for read in infile.fetch(region["chrom"]):
            outfile.write(read)


# In[4]: Main: Subset CRAM/BAM
def subset_cram(
    input_cram: str,
    output_file: str,
    region: Dict[str, Optional[int]],
    input_crai: Optional[str] = None,
    reference_fasta: Optional[str] = None,
    threads: int = 1
):
    """
    Subset a CRAM or BAM file to a specific region and save as BAM or CRAM.
    
    Args:
        input_cram: Path to input CRAM or BAM file
        output_file: Path to output file (.bam or .cram)
        region: Dictionary with 'chrom', 'start', 'end' keys
        input_crai: Optional CRAM index file path
        reference_fasta: Optional reference FASTA file path
        threads: Number of threads to use for compression/decompression
    """
    out_mode = determine_mode(output_file)
    open_kwargs = build_open_kwargs(reference_fasta, input_crai)
    infile, outfile = open_files(input_cram, output_file, out_mode, open_kwargs, threads)

    write_subset(infile, outfile, region)

    infile.close()
    outfile.close()

    if out_mode == "wb":  # Index BAM files only
        print(f"Indexing BAM file with {threads} threads...")
        pysam.index(output_file, "-@", str(threads))

    print(f"Subset file written to {output_file} using {threads} threads")
    return output_file

# In[5]: Helper: Get available CPU count
def get_available_cpus() -> int:
    """
    Get the number of available CPUs on the system.
    Returns the number of logical CPUs available.
    """
    try:
        return multiprocessing.cpu_count()
    except NotImplementedError:
        # Fallback if cpu_count() is not available
        return 1

# In[6]: Entry point for command line execution 
if __name__ == "__main__":
    # Get available CPUs for help text
    available_cpus = get_available_cpus()
    
    parser = argparse.ArgumentParser(
        description="Subset a CRAM or BAM file to a specific region and save as BAM or CRAM.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=f"""
Examples:
  %(prog)s input.cram output.bam chr6 28000000 34000000
  %(prog)s input.bam output.cram chr6 28000000 34000000 --threads 8
  %(prog)s input.cram output.bam chr6 28000000 34000000 --reference ref.fa --crai input.cram.crai

System Info:
  Available CPUs: {available_cpus}
        """
    )

    parser.add_argument("--input_file", required=True, help="Path to input CRAM or BAM file")
    parser.add_argument("--output_file", required=True, help="Path to output file (.bam or .cram)")
    parser.add_argument(
        "--chrom", type=int, choices=range(1, 24), required=True,
        help="Chromosome number (1-23). Will be converted to 'chrN'."
    )
    parser.add_argument("--start", type=int, default=None, help="Optional: Start position of region")
    parser.add_argument("--end", type=int, default=None, help="Optional: End position of region")
    parser.add_argument("--crai", default=None, help="Optional CRAM index file")
    parser.add_argument("--reference", default=None, help="Optional reference FASTA file")
    parser.add_argument("--threads", "-t", type=int, default=available_cpus,
                        help=f"Number of threads (default: {available_cpus}, auto-detected)")

    args = parser.parse_args()


    # Validate thread count
    if args.threads < 1:
        print(f"Error: Number of threads must be at least 1", file=sys.stderr)
        sys.exit(1)
    
    if args.threads > available_cpus:
        print(f"Warning: Requested {args.threads} threads, but only {available_cpus} CPUs available", 
              file=sys.stderr)

    # Format chrom as string for pysam
    chrom_str = f"chr{args.chrom}"

    region: Dict[str, Optional[int]] = {
        "chrom": chrom_str,
        "start": args.start,
        "end": args.end
    }


    print(f"Starting subset operation with {args.threads} threads...")
    if args.start and args.end:
        print(f"Region: {chrom_str}:{args.start}-{args.end}")
    else:
        print(f"Region: full {chrom_str}")

    subset_cram(
        input_cram=args.input_file,
        output_file=args.output_file,
        region=region,
        input_crai=args.crai,
        reference_fasta=args.reference,
        threads=args.threads
    )