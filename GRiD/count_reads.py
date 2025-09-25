# In[0]: Imports
import argparse
import pysam
from pathlib import Path
from GRiD.utils.utils import format_region, arg_check, argument_parser

# In[1]: Functions
def count_reads(bam_or_cram: str, region: str, reference: str | None = None) -> int:
    """
    Count reads overlapping the specified region.

    Parameters
    ----------
    bam_or_cram : str
        Path to BAM or CRAM file (must be indexed).
    region : str
        Genomic region string ("chr:start-end").
    reference : str | None
        Path to reference FASTA (required for CRAM unless reference is embedded).

    Returns
    -------
    int
        Number of reads overlapping the region.
    """
    mode = "rc" if bam_or_cram.endswith(".cram") else "rb"
    with pysam.AlignmentFile(bam_or_cram, mode, reference_filename=reference) as aln:
        return sum(1 for _ in aln.fetch(region=region))


# In[2]: Main
def main():
    parser = argument_parser()
    
    args = parser.parse_args()

    arg_check(args)

    if args.region is None:
        args.region = format_region(args)

    sample_id = Path(args.input).stem
    read_count = count_reads(args.input, args.region, reference=args.reference)

    with open(args.output, "w") as out:
        out.write("sample\tregion\tcount\n")
        out.write(f"{sample_id}\t{args.region}\t{read_count}\n")

# In[3]: Run
if __name__ == "__main__":
    main()