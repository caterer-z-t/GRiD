# In[0]: Imports
from pathlib import Path
from subset_cram import subset_bam_or_cram
from count_reads import count_reads
from mosdepth import run_mosdepth
from mosdepth_roi import extract_region_cov

# In[1]: Function to process a directory
def process_directory(input_dir, output_dir, reference, region, window=1000):
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True, parents=True)

    results = []

    for input_file in input_dir.glob("*.[bc]ram"):
        sample_id = input_file.stem
        print(f"\nProcessing {sample_id} ...")

        # subset
        subset_file = output_dir / f"{sample_id}.subset.bam"
        subset_bam_or_cram(str(input_file), region, str(subset_file), reference=reference)

        # count reads
        n_reads = count_reads(str(subset_file), region, reference=reference)

        # mosdepth
        mosdepth_prefix = output_dir / f"{sample_id}.mosdepth"
        run_mosdepth(subset_file, mosdepth_prefix, Path(reference), window=window)

        # extract region coverage
        regions_bed_gz = mosdepth_prefix.with_suffix(".regions.bed.gz")
        coverage = extract_region_cov(regions_bed_gz, region)

        results.append({
            "sample": sample_id,
            "reads": n_reads,
            "coverage": coverage
        })

    return results

# In[2]: Run example
if __name__ == "__main__":
    INPUT_DIR = "/path/to/cram_files"
    OUTPUT_DIR = "/path/to/output"
    REFERENCE = "/path/to/reference.fa"
    REGION = "chr6:160605062-160605662"

    all_results = process_directory(INPUT_DIR, OUTPUT_DIR, REFERENCE, REGION)
    for r in all_results:
        print(r)
