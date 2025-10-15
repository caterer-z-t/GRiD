# In[0]: Imports
import argparse
import pysam
from pathlib import Path
from GRiD.utils.utils import format_region, arg_check, argument_parser, ensure_index

# In[1]: Functions
def subset_bam_or_cram(
    bam_or_cram: str, region: str, output: str, reference: str | None = None
) -> None:
    """
    Subset a BAM/CRAM file to the specified region and save to a new file.

    Parameters
    ----------
    bam_or_cram : str
        Path to BAM or CRAM file (must be indexed).
    region : str
        Genomic region string ("chr:start-end").
    output : str
        Path to output BAM or CRAM file.
    reference : str | None
        Path to reference FASTA (required for CRAM unless reference is embedded).
    """
    mode_in = "rc" if bam_or_cram.endswith(".cram") else "rb"
    mode_out = "wc" if output.endswith(".cram") else "wb"

    with pysam.AlignmentFile(bam_or_cram, mode_in, reference_filename=reference) as aln_in, \
         pysam.AlignmentFile(output, mode_out, header=aln_in.header, reference_filename=reference) as aln_out:
        for read in aln_in.fetch(region=region):
            aln_out.write(read)


# In[2]: Main
def main():
    parser = argument_parser()
    args = parser.parse_args()
    arg_check(args)

    if args.region is None:
        args.region = format_region(args)

    print(f"Subsetting {args.input} for region {args.region}...")
    subset_bam_or_cram(
        args.input,
        args.region,
        args.output,
        reference=args.reference
    )

    print("Indexing output file...")
    # ensure output index is created
    ensure_index(args.output)

    print(f"Subset written to {args.output} for region {args.region}")


# In[3]: Run
if __name__ == "__main__":
    main()
