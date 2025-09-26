# In[0]: Imports
import gzip
from pathlib import Path
from GRiD.utils.utils import format_region, arg_check, argument_parser

# In[1]: Functions
def extract_region_cov(regions_bed_gz: Path, region: str) -> int:
    """
    Extract coverage for a given region from a mosdepth .regions.bed.gz file.

    Parameters
    ----------
    regions_bed_gz : Path
        Path to mosdepth .regions.bed.gz file.
    region : str
        Genomic region string ("chr:start-end").

    Returns
    -------
    int
        Coverage scaled by 100 and rounded (like int(100*$4+0.5) in bash).
    """
    chrom, coords = region.split(":")
    start, end = map(int, coords.split("-"))

    with gzip.open(regions_bed_gz, "rt") as f:
        for line in f:
            fields = line.strip().split("\t")
            if len(fields) < 4:
                continue
            c, s, e, cov = fields[:4]
            if c == chrom and int(s) == start and int(e) == end:
                return int(round(float(cov) * 100))

    raise ValueError(f"Region {region} not found in {regions_bed_gz}")


# In[2]: Main
def main():
    parser = argument_parser()
    args = parser.parse_args()
    arg_check(args)

    if args.region is None:
        args.region = format_region(args)

    # regions.bed.gz is expected to share prefix with mosdepth output
    regions_bed_gz = Path(args.output).with_suffix(".regions.bed.gz")

    if not regions_bed_gz.exists():
        raise FileNotFoundError(f"{regions_bed_gz} does not exist. Run mosdepth first.")

    sample_id = Path(args.input).stem
    coverage = extract_region_cov(regions_bed_gz, args.region)

    with open(args.output, "a") as out:
        out.write(f"{sample_id}\t{coverage}\n")


# In[3]: Run
if __name__ == "__main__":
    main()