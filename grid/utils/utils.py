# In[0]: Imports
import os
import pysam    
from pathlib import Path
from threading import Lock
import gzip
import glob
from rich.progress import Progress, SpinnerColumn, TextColumn, BarColumn, TaskProgressColumn
from contextlib import contextmanager

# In[0.1]: Utility functions
def log(console, msg, style=None):
    if console:
        if style:
            console.print(msg, style=style)
        else:
            console.print(msg)
    else:
        print(msg)

@contextmanager
def progress_bar(console=None, total=1, description="Working"):
    """
    Reusable progress bar context manager.
    
    Usage:
        with progress_bar(console, total=len(samples), description="Processing samples") as (progress, task):
            for sample in samples:
                progress.update(task, description=f"Processing {sample}")
                # do work
                progress.advance(task)
    """
    with Progress(
        SpinnerColumn(spinner_name="dots", style="info"),
        TextColumn("[progress.description]{task.description}", style="highlight"),
        BarColumn(complete_style="success", finished_style="success"),
        TaskProgressColumn(),
        console=console
    ) as progress:
        task = progress.add_task(description, total=total)
        yield progress, task


def find_file(directory_loc, sample, expected_type=None):
    """Return file_path or None."""
    if expected_type:
        pattern = os.path.join(directory_loc, f"*{sample}*.{expected_type}")
        matches = glob.glob(pattern)
        if matches:
            return matches[0]
    return None

def has_index(file_path, file_type):
    """Check if appropriate index exists."""

    allowed_file_types = {
        "CRAM": "crai",
        "BAM": "bai"
    }
    if file_type.upper() not in allowed_file_types.keys():
        return False
    
    if file_type.upper() == "CRAM":
        return ( 
            os.path.exists(file_path + ".crai") or 
            os.path.exists(file_path.replace(".cram", ".crai")) 
        )
    
    if file_type.upper() == "BAM":
        return (
            os.path.exists(file_path + ".bai") or
            os.path.exists(file_path.replace(".bam", ".bai"))
        )
    
    return False

def get_samples(samples_file):
    with open(samples_file) as f:
        return [line.strip() for line in f if line.strip()]

def get_flags(config, parameter):
    return [flag for flag in config.get(parameter, {}).get('flags', []) if flag is not None]

def create_index_for_file(file_path, file_type, reference_genome):
    if file_type.upper() == "CRAM":
        pysam.index(file_path, file_path + ".crai", reference_filename=reference_genome)
    elif file_type.upper() == "BAM":
        pysam.index(file_path, file_path + ".bai", reference_filename=reference_genome)

def setup_output_file(output_file: str, chrom: str, start: int, end: int) -> Path:
    """
    Create output file with header.
    
    Args:
        output_file: Path to output TSV file
        chrom: Chromosome name
        start: Start position
        end: End position
    
    Returns:
        Path object for the output file
    """
    output_path = Path(output_file).expanduser()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    with open(output_path, "w") as f:
        f.write(f"Sample\t{chrom}:{start}-{end}\n")
    
    return output_path
# In[1]: Helper function to check for BAM/CRAM files and their indexes
def check_index(
        config,
        console=None
    ):
    """Check that all samples have BAM/CRAM files and corresponding indexes."""
    
    try:
        file_type = config.get('file_type')
        directory_loc = config['directory_loc']
        samples_file = config['samples_file']
        output_dir = config.get('output_dir', '.')
    except Exception as e:
        log(console, f"Config error: {e}", style="danger")
        return

    results = {
        "missing_file": [],
        "missing_index": [],
        "has_index": []
    }

    # ---- main ----
    samples = get_samples(samples_file)

    with progress_bar(console, total=len(samples), description="Checking indexes") as (progress, task):
        for sample in samples:
            progress.update(task, description=f"Checking {sample}", style="info")
            file_path = find_file(directory_loc, sample, file_type)

            if not file_path:
                results["missing_file"].append(sample)
                progress.advance(task)

            if has_index(file_path, file_type):
                results["has_index"].append(sample)
                progress.advance(task)

            else:
                results["missing_index"].append(sample)
                progress.advance(task)

    if config['index'].get('output_file_prefix'):
        output_file = f"{output_dir}/{config['index']['output_file_prefix']}.{config.get('output_file_type', 'tsv')}"
        with open(output_file, "w") as f:
            f.write("Sample\tStatus\n")
            for sample in results["has_index"]:
                f.write(f"{sample}\tHas index\n")
            for sample in results["missing_file"]:
                f.write(f"{sample}\tMissing file\n")
            for sample in results["missing_index"]:
                f.write(f"{sample}\tMissing index\n")

        log(console, f"Index check results written to {output_file}", style="success")

# In[2]: Helper function to create index for BAM/CRAM files
def create_index(
        config, 
        console=None
    ):
    """Create index for BAM/CRAM files if not already present."""
    try:
        file_type = config.get('file_type')
        directory_loc = config['directory_loc']
        samples_file = config['samples_file']
        reference_genome = config['reference_genome']
    except Exception as e:
        log(console, f"Config error: {e}", style="danger")
        return
    
    results = {
        "missing_file": [],
        "missing_index": [],
        "has_index": []
    }
        
    # ---- main ----
    samples = get_samples(samples_file)

    with progress_bar(console, total=len(samples), description="Creating index") as (progress, task):
        for sample in samples:
            progress.update(task, description=f"Processing {sample}", style="info")
            file_path = find_file(directory_loc, sample, file_type)

            if not file_path:
                results["missing_file"].append(sample)
                progress.advance(task)
                
                continue

            if has_index(file_path, file_type):
                results["has_index"].append(sample)
                progress.advance(task)
                continue

            try:
                create_index_for_file(file_path, file_type, reference_genome)
            except Exception as e:
                log(console, f"Failed to create index for {sample}: {e}", style="danger")
                results["missing_index"].append(sample)

            progress.advance(task)

    if config['index'].get('output_file_prefix') and (results["missing_file"] or results["missing_index"]):
        output_file = f"{config.get('output_dir', '.')}/{config['index']['output_file_prefix']}.err"
        with open(output_file, "w") as f:
            f.write("Sample\tStatus\n")
            for sample in results["has_index"]:
                f.write(f"{sample}\tHas index\n")
            for sample in results["missing_file"]:
                f.write(f"{sample}\tMissing file\n")
            for sample in results["missing_index"]:
                f.write(f"{sample}\tFailed to create index\n")

        log(console, f"Index creation results written to {output_file}", style="success")

def create_region_string(chrom: str, start: int, end: int, console=None) -> str:
    """
    Create a genomic region string in 'chr:start-end' format.

    Args:
        chrom: Chromosome name (e.g., 'chr6')
        start: Start position (e.g., 160000000)
        end: End position (e.g., 160100000)

    Returns:
        Genomic region string.
    """
    
    if chrom and start is not None and end is not None:
        region = f"{chrom}:{start}-{end}"
    else:
        log(console, "You must provide either a full region or chromosome, start, and end separately.", style="danger")
        raise ValueError("Invalid region parameters")
    return region


def open_maybe_gz(path, mode="rt"):
    if str(path).endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)