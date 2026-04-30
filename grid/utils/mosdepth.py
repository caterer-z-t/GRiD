# In[0]: Imports
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from threading import Lock
from functools import partial
import shutil
import subprocess
import time
from typing import Tuple, List
import gzip

from .utils import log, get_samples, setup_output_file, find_file, progress_bar


# In[1]: Main Mosdepth
def compute_mosdepth(config, console=None):
    try:
        directory_loc = config["directory_loc"]
        samples_file = config["samples_file"]
        samples = get_samples(samples_file)
        chrom = config.get("chrom", None)
        start = config.get("start_bp", None)
        end = config.get("end_bp", None)

        output_file_prefix = config.get("mosdepth", {}).get("output_file_prefix", None)
        output_file_type = config.get("output_file_type", "tsv")
        output_dir = config.get("output_dir", ".")
        output_file = Path(f"{output_dir}/{output_file_prefix}.{output_file_type}")

        threads = config.get("threads", 1)
        ref = config.get("reference_genome", {})
        region_name = config.get("mosdepth", {}).get("region_name", None)
        by = config.get("mosdepth", {}).get("bin_size", 1000)
        fast_mode = config.get("mosdepth", {}).get("mode", "fast").lower()
        work_dir = config.get("mosdepth", {}).get("work_dir", None)
        remove_intermediate = config.get("mosdepth", {}).get("remove_intermediate", False)
    except Exception as e:
        log(console, f"Config error: {e}", style="danger")
        return

    # Check mosdepth availability
    try:
        check_mosdepth_available()
    except RuntimeError as e:
        log(console, f"✗ {str(e)}", style="danger")
        return

    work_path = Path(work_dir).expanduser()
    work_path.mkdir(parents=True, exist_ok=True)

    output_file = setup_output_file(output_file, chrom, start, end)

    files = {
        sample: result
        for sample in samples
        if (result := find_file(directory_loc, sample, config.get("file_type")))[0] is not None
    }

    # Create lock for thread-safe file writing
    write_lock = Lock()

    # Create partial function with fixed parameters
    process_func = partial(
        run_mosdepth_single_cram,
        ref_fasta=ref,
        work_dir=work_path,
        chrom=chrom,
        start=start,
        end=end,
        region_name=region_name,
        by=by,
        fast_mode=fast_mode,
        threads=threads,
        console=console,
    )

    # Process CRAMs with threading and progress tracking
    failed = []

    with progress_bar(console, total=len(files), description="Running mosdepth") as (
        progress,
        task,
    ):
        with ThreadPoolExecutor(max_workers=threads) as executor:
            # Submit all jobs
            future_to_file = {
                executor.submit(process_func, file): sample for sample, file in files.items()
            }

            # Collect and write results as they complete
            for future in as_completed(future_to_file):
                sample = future_to_file[future]
                coverage = future.result()

                if coverage != "Error":
                    write_coverage_result(output_file, sample, coverage, write_lock)
                else:
                    log(console, f"✗ {sample} failed", style="danger")
                    failed.append(sample)

                progress.update(task, advance=1)

    log(console, f"Mosdepth coverage results written to {output_file}", style="success")
    if remove_intermediate:
        remove_intermediate_files(work_path, console)


# In[2]: Helper Functions
def check_mosdepth_available():
    """
    Check if mosdepth executable is available in PATH.

    Args:
        None

    Returns:
        None

    Raises:
        RuntimeError: If mosdepth is not found.
    """
    if shutil.which("mosdepth") is None:
        raise RuntimeError("mosdepth not found in PATH. Please install mosdepth before running.")


def write_coverage_result(
    output_file: Path, sample_name: str, coverage: int, write_lock: Lock
) -> None:
    """
    Write coverage result to output file in a thread-safe manner.

    Args:
        output_file: Path to output TSV file
        sample_name: Sample name
        coverage: Coverage value
        write_lock: Threading lock for safe file writing
    """
    with write_lock:
        with open(output_file, "a", newline="") as f:
            f.write(f"{sample_name}\t{coverage}\n")


def run_mosdepth_single_cram(
    cram_path: str,
    ref_fasta: str,
    work_dir: Path,
    chrom: str,
    start: int,
    end: int,
    region_name: str,
    by: int,
    fast_mode: bool,
    threads: int = 1,
    console=None,
) -> Tuple[str, int | str]:
    """
    Run mosdepth on a single CRAM file and compute coverage for region.

    Args:
        cram_path: Path to CRAM file
        ref_fasta: Path to reference genome FASTA
        work_dir: Working directory for intermediate files
        chrom: Chromosome name
        start: Start position
        end: End position
        by: Bin size for mosdepth
        fast_mode: Whether to use fast mode
        threads: Number of threads to use
    Returns:
        Tuple of (sample_name, coverage) where coverage is int or "Error"
    """
    cram = Path(cram_path)
    sample_name = cram.stem

    try:
        # Build and run mosdepth command
        out_prefix = work_dir / f"{sample_name}_{region_name}"
        cmd = build_mosdepth_command(str(cram), ref_fasta, out_prefix, by, fast_mode, threads)
        subprocess.run(cmd, check=True, capture_output=True, text=True)

        # Wait for output file
        regions_file = wait_for_mosdepth_output(work_dir, sample_name, console)

        # Compute coverage
        coverage = compute_region_coverage(regions_file, chrom, start, end)

        return coverage

    except Exception as e:
        return "Error"


def build_mosdepth_command(
    cram_path: str, ref_fasta: str, output_prefix: Path, by: int, fast_mode: bool, threads: int = 1
) -> List[str]:
    """
    Build mosdepth command with appropriate flags.

    Args:
        cram_path: Path to CRAM file
        ref_fasta: Path to reference genome FASTA
        output_prefix: Output prefix for mosdepth
        by: Bin size for coverage calculation
        fast_mode: Whether to use fast mode
        threads: Number of threads to use
    Returns:
        List of command arguments
    """
    cmd = [
        "mosdepth",
        "-n",
        "--by",
        str(by),
        "-f",
        str(ref_fasta),
        str(output_prefix),
        str(cram_path),
        "-t",
        str(threads),
    ]

    if fast_mode:
        cmd.insert(1, "--fast-mode")

    return cmd


def wait_for_mosdepth_output(
    work_dir: Path, sample_name: str, console=None, max_attempts: int = 3, sleep_seconds: int = 2
) -> Path:
    """
    Wait for mosdepth output file to appear.

    Args:
        work_dir: Working directory where mosdepth writes output
        sample_name: Sample name to search for
        max_attempts: Maximum number of attempts to find file
        sleep_seconds: Seconds to wait between attempts

    Returns:
        Path to regions.bed.gz file

    Raises:
        FileNotFoundError: If file is not found after max attempts
    """
    for attempt in range(max_attempts):
        matches = list(work_dir.glob(f"{sample_name}*regions.bed.gz"))
        if matches:
            return matches[0]

        log(
            console,
            f"Waiting for mosdepth output ({attempt + 1}/{max_attempts})...",
            style="warning",
        )
        time.sleep(sleep_seconds)

    raise FileNotFoundError(
        f"Expected mosdepth output matching pattern '{sample_name}*regions.bed.gz' "
        f"not found in {work_dir} for sample {sample_name}."
    )


def compute_region_coverage(regions_file: Path, chrom: str, start: int, end: int) -> int:
    """
    Compute average coverage for a genomic region from mosdepth output.

    Args:
        regions_file: Path to mosdepth regions.bed.gz file
        chrom: Chromosome name
        start: Start position
        end: End position

    Returns:
        Average coverage as integer (rounded and scaled by 100)
    """
    region_cov = 0.0
    covered_bp = 0

    with gzip.open(regions_file, "rt") as f:
        for line in f:
            r_chr, r_start, r_end, mean_cov = line.strip().split("\t")[:4]
            r_start, r_end = int(r_start), int(r_end)
            mean_cov = float(mean_cov)

            if r_chr != chrom:
                continue

            overlap_start = max(start, r_start)
            overlap_end = min(end, r_end)
            overlap = overlap_end - overlap_start

            if overlap > 0:
                region_cov += mean_cov * overlap
                covered_bp += overlap

    return int(round(100 * (region_cov / covered_bp))) if covered_bp > 0 else 0


def remove_intermediate_files(
    work_dir: Path, console=None, include_region_bed_gz: bool = False
) -> None:
    """
    Remove intermediate mosdepth files from the working directory.

    Args:
        work_dir: Path to working directory
        console: Optional console for logging
        include_region_bed_gz: Whether to include regions.bed.gz files in removal

    Returns:
        None
    """
    files_to_remove = ["mosdepth.global.dist.txt", "mosdepth.region.dist.txt", "regions.bed.gz.csi"]

    work_dir = Path(work_dir)

    if include_region_bed_gz:
        files_to_remove.append("regions.bed.gz")

    for file in work_dir.glob("*"):
        if any(file.name.endswith(suffix) for suffix in files_to_remove):
            try:
                file.unlink()
            except Exception as e:
                log(console, f"Failed to remove intermediate file {file}: {e}", style="warning")
