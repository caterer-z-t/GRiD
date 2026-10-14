# In[0]: Imports
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from functools import partial
import shutil
import subprocess
import gzip

from .utils import log, get_samples, find_file, progress_bar


# In[1]: Main Mosdepth
def compute_mosdepth(config, console=None):
    try:
        directory_loc = config["directory_loc"]
        samples_file = config["samples_file"]
        samples = get_samples(samples_file)

        output_file_prefix = config.get("mosdepth", {}).get("output_file_prefix", "mosdepth_bins")
        output_file_type = config.get("output_file_type", "txt")
        output_dir = Path(config.get("output_dir", "."))
        threads = config.get("threads", 1)
        ref = config.get("reference_genome")
        by = config.get("mosdepth", {}).get("bin_size", 1000)
        mosdepth_dir = Path(config.get("mosdepth", {}).get("work_dir", output_dir / "mosdepth"))
        remove_intermediate = config.get("mosdepth", {}).get("remove_intermediate", True)

        # Example bed: use configured path if provided, otherwise default to mosdepth_dir
        example_bed_path = config.get("mosdepth", {}).get("example_bed", None)
        if example_bed_path is not None:
            example_bed = Path(example_bed_path)
        else:
            example_bed = mosdepth_dir / "example.regions.bed.gz"

        combined_file = output_dir / f"{output_file_prefix}.{output_file_type}"
        combined_gz = Path(str(combined_file) + ".gz")

    except Exception as e:
        log(console, f"Config error: {e}", style="danger")
        return

    # Check mosdepth availability
    if shutil.which("mosdepth") is None:
        log(console, "mosdepth not found in PATH. Please install mosdepth.", style="danger")
        return

    mosdepth_dir.mkdir(parents=True, exist_ok=True)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Remove existing combined/gz files if present so we start fresh
    combined_file.unlink(missing_ok=True)
    combined_gz.unlink(missing_ok=True)

    files = {
        sample: result
        for sample in samples
        if (result := find_file(directory_loc, sample, config.get("file_type"))) is not None
    }

    if not files:
        log(console, "No CRAM files found — check directory_loc and samples_file.", style="danger")
        return

    # Process each CRAM: run mosdepth, parse regions.bed.gz, write per-sample file
    process_func = partial(
        run_mosdepth_single_cram,
        ref_fasta=ref,
        mosdepth_dir=mosdepth_dir,
        by=by,
        console=console,
    )

    failed = []
    with progress_bar(console, total=len(files), description="Running mosdepth") as (
        progress,
        task,
    ):
        with ThreadPoolExecutor(max_workers=threads) as executor:
            future_to_sample = {
                executor.submit(process_func, cram_path): sample
                for sample, cram_path in files.items()
            }

            for future in as_completed(future_to_sample):
                sample = future_to_sample[future]
                success = future.result()
                if not success:
                    failed.append(sample)
                    log(console, f"✗ {sample} failed", style="danger")
                progress.advance(task)

    if failed:
        log(console, f"{len(failed)} samples failed: {failed}", style="warning")

    # Save example bed from first successful sample if not already present
    if not example_bed.exists():
        per_sample_beds = sorted(mosdepth_dir.glob("*.regions.bed.gz"))
        if per_sample_beds:
            shutil.copy(per_sample_beds[0], example_bed)
            log(console, f"Saved example regions bed: {example_bed}", style="info")
        else:
            log(
                console,
                "WARNING: No regions.bed.gz files found to save as example.",
                style="warning",
            )
    else:
        log(
            console,
            f"Example regions bed already exists — skipping: {example_bed}",
            style="info",
        )

    # Verify all expected per-sample files exist before merging
    per_sample_files = sorted(mosdepth_dir.glob("*.per_sample.txt"))
    n_expected = len(files)
    n_found = len(per_sample_files)
    if n_found < n_expected:
        log(
            console,
            f"WARNING: Expected {n_expected} per-sample files, found {n_found}. "
            f"{n_expected - n_found} samples may have failed.",
            style="warning",
        )
    log(console, f"Merging {n_found} per-sample files into combined output...", style="info")

    # Merge all per-sample files into combined file sequentially (no threads)
    with open(combined_file, "w") as out_f:
        for sample_file in per_sample_files:
            out_f.write(sample_file.read_text())

    # Verify merged file looks correct before gzipping
    n_lines = sum(1 for _ in open(combined_file))
    if n_lines != n_found:
        log(
            console,
            f"ERROR: Merged file has {n_lines} lines but expected {n_found} — aborting gzip.",
            style="danger",
        )
        return
    log(console, f"Merged file has {n_lines} lines — looks correct.", style="info")

    # Gzip using system gzip (more reliable than Python's gzip for large files)
    log(console, "Gzipping combined file...", style="info")
    try:
        subprocess.run(
            ["gzip", "-f", str(combined_file)],
            check=True,
            capture_output=True,
            text=True,
        )
    except subprocess.CalledProcessError as e:
        log(console, f"gzip failed: {e.stderr}", style="danger")
        return

    # Verify gz file exists and is not empty before declaring success
    if not combined_gz.exists() or combined_gz.stat().st_size == 0:
        log(console, f"ERROR: gzip output missing or empty: {combined_gz}", style="danger")
        return

    # Spot-check gz integrity
    try:
        result = subprocess.run(
            ["gzip", "-t", str(combined_gz)],
            check=True,
            capture_output=True,
            text=True,
        )
    except subprocess.CalledProcessError:
        log(console, f"ERROR: gzip integrity check failed for {combined_gz}", style="danger")
        return

    # Remove intermediate files if requested — only after gz is confirmed good
    if remove_intermediate:
        remove_intermediate_files(mosdepth_dir, console)

    log(console, f"Done. Combined output: {combined_gz}", style="success")
    log(console, f"Example regions bed: {example_bed}", style="success")


# In[2]: Run mosdepth on a single CRAM
def run_mosdepth_single_cram(
    cram_path: str,
    ref_fasta: str,
    mosdepth_dir: Path,
    by: int = 1000,
    console=None,
) -> bool:
    """
    Run mosdepth genome-wide on one CRAM, parse the regions.bed.gz output
    exactly as Mukamel's shell script does (strip chr, keep autosomes 1-22,
    scale depth by 100 and round to int), and write one per-sample row file.

    Args:
        cram_path: Path to CRAM file
        ref_fasta: Path to reference FASTA
        mosdepth_dir: Directory for mosdepth intermediate output
        by: Bin size (default 1000 = 1kb)
        console: Optional rich console for logging

    Returns:
        True on success, False on failure
    """
    cram = Path(cram_path)
    sample_id = cram.stem
    out_prefix = mosdepth_dir / sample_id
    per_sample_out = mosdepth_dir / f"{sample_id}.per_sample.txt"

    # Skip if already processed
    if per_sample_out.exists():
        return True

    try:
        # Run mosdepth genome-wide, no region subsetting, fast mode, no per-base output
        cmd = [
            "mosdepth",
            "-n",
            "--fast-mode",
            "--by",
            str(by),
            "-f",
            str(ref_fasta),
            str(out_prefix),
            str(cram),
        ]
        subprocess.run(cmd, check=True, capture_output=True, text=True)

        regions_file = mosdepth_dir / f"{sample_id}.regions.bed.gz"

        if not regions_file.exists():
            log(console, f"regions.bed.gz not found for {sample_id}", style="danger")
            return False

        # Parse regions.bed.gz exactly as the awk in the shell script:
        # strip chr prefix, keep autosomes 1-22, scale by 100 and round
        row_values = [sample_id]
        with gzip.open(regions_file, "rt") as f:
            for line in f:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 4:
                    continue
                c = parts[0]
                if c.startswith("chr"):
                    c = c[3:]
                # Keep only autosomes 1-22 (length<=2 and numeric 1-22)
                if len(c) <= 2:
                    try:
                        c_int = int(c)
                    except ValueError:
                        continue
                    if 1 <= c_int <= 22:
                        depth = float(parts[3])
                        row_values.append(str(int(100 * depth + 0.5)))

        # Write sample row: "SAMPLE_ID\tval1\tval2\t...\n"
        with open(per_sample_out, "w") as f:
            f.write("\t".join(row_values) + "\n")

        return True

    except Exception as e:
        log(console, f"Error processing {sample_id}: {e}", style="danger")
        return False


# In[3]: Remove intermediate files
def remove_intermediate_files(mosdepth_dir: Path, console=None) -> None:
    """
    Remove intermediate mosdepth files, keeping only the example bed.

    Removes:

    - ``*.mosdepth.global.dist.txt``
    - ``*.mosdepth.region.dist.txt``
    - ``*.regions.bed.gz``
    - ``*.regions.bed.gz.csi``
    - ``*.per_sample.txt``
    """

    suffixes = [
        ".mosdepth.global.dist.txt",
        ".mosdepth.region.dist.txt",
        ".regions.bed.gz",
        ".regions.bed.gz.csi",
        ".per_sample.txt",
    ]

    for f in mosdepth_dir.glob("*"):
        if f.name == "example.regions.bed.gz":
            continue
        if any(f.name.endswith(s) for s in suffixes):
            try:
                f.unlink()
            except Exception as e:
                log(console, f"Failed to remove {f}: {e}", style="warning")
