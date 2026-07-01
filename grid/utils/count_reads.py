# In[0]: Imports
from __future__ import annotations

from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from threading import Lock
from functools import partial
import pysam

from .utils import log, get_samples, setup_output_file, find_file, progress_bar


# In[1]: Main count reads
def count_reads(config, console=None):
    try:
        directory_loc = config["directory_loc"]
        samples_file = config["samples_file"]
        samples = get_samples(samples_file)
        chrom = config.get("chrom", None)
        start = config.get("start_bp", None)
        end = config.get("end_bp", None)
        flags = config.get("count_reads", {}).get("flags", [])
        threads = config.get("threads", 1)
        min_mapq = config.get("min_mapq", 1)

        output_file_prefix = config.get("count_reads", {}).get("output_file_prefix", None)
        output_file_type = config.get("output_file_type", "tsv")
        output_dir = config.get("output_dir", ".")
        output_file = Path(f"{output_dir}/{output_file_prefix}.{output_file_type}")

        ref = config.get("reference_genome")
    except Exception as e:
        log(console, f"Config error: {e}", style="danger")
        return

    # Setup output file with header
    output_path = setup_output_file(output_file, chrom, start, end)

    files = {
        sample: result
        for sample in samples
        if (result := find_file(directory_loc, sample, config.get("file_type")))[0] is not None
    }

    # Create lock for thread-safe file writing
    write_lock = Lock()

    # Create a partial function with fixed parameters
    process_func = partial(
        process_single_cram,
        ref_fasta=ref,
        chrom=chrom,
        start=start,
        end=end,
        proper_flags=flags,
        min_mapq=min_mapq,
        console=console,
    )

    # Process CRAMs with threading and progress tracking
    with progress_bar(console, total=len(files), description="Counting reads") as (progress, task):
        with ThreadPoolExecutor(max_workers=threads) as executor:
            # Submit all jobs
            future_to_sample = {
                executor.submit(process_func, file): sample for sample, file in files.items()
            }

            # Collect and write results as they complete
            for future in as_completed(future_to_sample):
                sample = future_to_sample[future]
                count = future.result()

                # Write result to file (thread-safe)
                write_read_results(output_path, sample, count, write_lock)

                # Update progress bar
                progress.advance(task)

    log(console, f"Read counting completed. " f"Results written to {output_path}", style="success")


def count_reads_in_region(
    cram_file: str,
    ref_fasta: str,
    chrom: str,
    start: int,
    end: int,
    proper_flags: set[int],
    min_mapq: int = 1,
) -> int:
    """
    Count properly paired reads in a genomic region using pysam.
    """
    count = 0
    with pysam.AlignmentFile(cram_file, "rc", reference_filename=ref_fasta) as bam:
        for read in bam.fetch(chrom, start, end):
            if (
                read.flag in proper_flags  # exact flag match like C++
                and read.mapq >= min_mapq
                and read.reference_id == read.next_reference_id
                and not read.is_duplicate  # exclude duplicates
                and not read.is_secondary  # exclude secondary alignments
                and read.reference_start
                >= start  # start position inside region (matchs Hujoel C++ bin logic)
                and read.reference_start < end
            ):
                count += 1
    return count


def process_single_cram(
    path: str,
    ref_fasta: str,
    chrom: str,
    start: int,
    end: int,
    proper_flags: set[int],
    min_mapq: int = 1,
    console=None,
) -> int | str:
    """
    Process a single sequence alignment file and count reads in region.

    Args:
        path: Path to sequence alignment file
        ref_fasta: Path to reference genome FASTA
        chrom: Chromosome name
        start: Start position
        end: End position
        proper_flags: Set of SAM flags to count

    Returns:
        int or str: Read count or "Error"
    """
    basename = Path(path).name

    try:
        count = count_reads_in_region(path, ref_fasta, chrom, start, end, proper_flags, min_mapq)
        return count
    except Exception as e:
        log(console, f"Failed to count reads for {basename}: {e}", style="danger")
        return "Error"


# In[3]: Helper function to write results to file in a thread-safe way
def write_read_results(
    output_file: Path, basename: str, count: int | str, write_lock: Lock
) -> None:
    """
    Write a single result to the output file in a thread-safe manner.

    Args:
        output_file: Path to output file
        basename: Sample basename
        count: Read count or "Error"
        write_lock: Threading lock for safe file writing
    """
    with write_lock:
        with open(output_file, "a") as f:
            f.write(f"{basename}\t{count}\n")
