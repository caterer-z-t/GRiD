# In[0]: Imports
from pathlib import Path
from GRiD.utils.utils import format_region, arg_check, argument_parser
from subset_cram import subset_bam_or_cram
from count_reads import count_reads
from mosdepth import run_mosdepth
from mosdepth_roi import extract_region_cov

# In[1]: Main function
def process_file(input_file, output_dir, reference, region=None, window=1000):
    input_path = Path(input_file)
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True, parents=True)
    sample_id = input_path.stem

    if region is None:
        raise ValueError("Region must be provided for single-file processing.")

    # Step 1: subset
    subset_file = output_dir / f"{sample_id}.subset.bam"
    subset_bam_or_cram(str(input_path), region, str(subset_file), reference=reference)
    print(f"Subset done: {subset_file}")

    # Step 2: count reads
    n_reads = count_reads(str(subset_file), region, reference=reference)
    print(f"Read count for {sample_id}: {n_reads}")

    # Step 3: mosdepth
    mosdepth_prefix = output_dir / f"{sample_id}.mosdepth"
    run_mosdepth(subset_file, mosdepth_prefix, Path(reference), window=window)
    print(f"Mosdepth done: {mosdepth_prefix}")

    # Step 4: extract region coverage
    regions_bed_gz = mosdepth_prefix.with_suffix(".regions.bed.gz")
    coverage = extract_region_cov(regions_bed_gz, region)
    print(f"Region coverage: {coverage}")

    return {
        "sample": sample_id,
        "reads": n_reads,
        "coverage": coverage
    }

# In[2]: CLI
def main():
    parser = argument_parser()
    args = parser.parse_args()
    arg_check(args)

    if args.region is None:
        args.region = format_region(args)

    results = process_file(
        args.input,
        args.output,
        reference=args.reference,
        region=args.region
    )
    print(results)

# In[3]: Run
if __name__ == "__main__":
    main()
