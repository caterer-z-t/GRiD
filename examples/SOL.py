#!/usr/bin/env python3
"""
SOL Pipeline (Python version)
=============================

Usage:
    python sol_pipeline.py        # runs steps 1–9
    python sol_pipeline.py 5      # runs only step 5
    python sol_pipeline.py 3-5    # runs steps 3–5

This script mirrors the functionality of the original SLURM script,
but calls GRiD functions directly from Python.
"""

import sys
import re
from pathlib import Path
from rich.console import Console
from grid.cli import (
    crai,
    count_reads,
    mosdepth,
    normalize_mosdepth,
    find_neighbors,
    extract_reference,
    lpa_realign,
    compute_dipcn,
    estimate_kiv,
)

console = Console()

# ============================================================
# Step Range Parsing
# ============================================================
def parse_range(arg: str):
    """Parse step range input (e.g., '3-5', '5', or '')."""
    if not arg:
        return 1, 9
    match = re.match(r"^(\d+)-(\d+)$", arg)
    if match:
        return int(match.group(1)), int(match.group(2))
    elif re.match(r"^\d+$", arg):
        step = int(arg)
        return step, step
    else:
        console.print(f"[red]❌ Invalid range argument: {arg}[/red]")
        sys.exit(1)


# ============================================================
# Environment Setup
# ============================================================
WORK = Path("/work/users/c/a/catererz/grid")
GRID = Path("/nas/longleaf/home/catererz/epi/GRiD")
SOL = Path("/proj/epi/Genetic_Data_Center/SOL")
CRAM = SOL / "cram"
FILES = GRID / "files"

WORK.mkdir(parents=True, exist_ok=True)

REGIONS_FILE = FILES / "734_possible_coding_vntr_regions.IBD2R_gt_0.25.uniq.txt"

# Parse LPA region
with REGIONS_FILE.open() as f:
    for line in f:
        cols = line.split()
        if len(cols) >= 7 and cols[6] == "LPA":
            CHR, START, END = cols[0], int(cols[1]), int(cols[2])
            break


# ============================================================
# Helper: should we run this step?
# ============================================================
def should_run(step, start, end):
    return start <= step <= end


# ============================================================
# Main pipeline
# ============================================================
def main():
    arg = sys.argv[1] if len(sys.argv) > 1 else ""
    start_step, end_step = parse_range(arg)
    console.rule(f"Running SOL Pipeline Steps {start_step}–{end_step}")

    # STEP 1 — Ensure CRAI Indexes
    if should_run(1, start_step, end_step):
        console.rule("[bold cyan]Step 1: Ensuring CRAM Indexes[/bold cyan]")
        for cram_file in CRAM.glob("*.cram"):
            crai.callback(cram=cram_file, reference=str(FILES / "hg38.fa"))

    # STEP 2 — Count Reads
    if should_run(2, start_step, end_step):
        console.rule("[bold cyan]Step 2: Counting Reads[/bold cyan]")
        count_reads.callback(
            cram_dir=str(CRAM),
            output_file="count-reads.tsv",
            ref_fasta=str(FILES / "hg38.fa"),
            chrom=f"chr{CHR}",
            start=START,
            end=END,
            config=str(GRID / "grid/config.yaml"),
            threads=4,
        )

    # STEP 3 — Mosdepth
    if should_run(3, start_step, end_step):
        console.rule("[bold cyan]Step 3: Running Mosdepth[/bold cyan]")
        mosdepth.callback(
            cram_dir=str(CRAM),
            output_file="mosdepth.tsv",
            ref_fasta=str(FILES / "hg38.fa"),
            chrom=f"chr{CHR}",
            start=START,
            end=END,
            region_name="LPA",
            work_dir=str(WORK / "mosdepth"),
            by=1000,
            fast=True,
            threads=4,
        )

    # STEP 4 — Normalize Mosdepth
    if should_run(4, start_step, end_step):
        console.rule("[bold cyan]Step 4: Normalizing Mosdepth[/bold cyan]")
        normalize_mosdepth.callback(
            mosdepth_dir=str(WORK / "mosdepth"),
            output_file="normalized_mosdepth.txt.gz",
            repeat_mask=str(FILES / "repeat_mask_list.hg38.ucsc_bed"),
            chrom=f"chr{CHR}",
            start=START,
            end=END,
            min_depth=20,
            max_depth=100,
            top_frac=0.1,
            threads=4,
        )

    # STEP 5 — Find Neighbors
    if should_run(5, start_step, end_step):
        console.rule("[bold cyan]Step 5: Finding Neighbors[/bold cyan]")
        find_neighbors.callback(
            input_file="normalized_mosdepth.txt.gz",
            output_file="neighbors.txt",
            zmax=2.0,
            n_neighbors=500,
            sigma2_max=1000.0,
        )

    # STEP 6 — Extract Reference
    if should_run(6, start_step, end_step):
        console.rule("[bold cyan]Step 6: Extracting LPA Reference[/bold cyan]")
        extract_reference.callback(
            reference_fa=str(FILES / "hg38.fa"),
            bed_file=str(FILES / "ref_hg38.bed"),
            output_dir=str(WORK),
            output_prefix="lpa_reference",
        )

    # STEP 7 — Realign Reads
    if should_run(7, start_step, end_step):
        console.rule("[bold cyan]Step 7: Realigning Reads to LPA[/bold cyan]")
        lpa_realign.callback(
            cram_dir=str(CRAM),
            output_file="lpa_realigned.tsv",
            ref_fasta=str(FILES / "hg38.fa"),
            lpa_ref_fasta="lpa_reference.fasta",
            positions_file=str(FILES / "hardcoded_positions.txt"),
            genome_build="hg38",
            chrom=f"chr{CHR}",
            start=START,
            end=END,
            threads=4,
        )

    # STEP 8 — Compute DIPCN
    if should_run(8, start_step, end_step):
        console.rule("[bold cyan]Step 8: Computing LPA DIPCN[/bold cyan]")
        compute_dipcn.callback(
            count_file="lpa_realigned.tsv",
            neighbor_file="neighbors.txt",
            output_prefix="lpa_dipcn",
            n_neighbors=500,
        )

    # STEP 9 — Compute KIV2 CN
    if should_run(9, start_step, end_step):
        console.rule("[bold cyan]Step 9: Computing LPA KIV2 CN[/bold cyan]")
        estimate_kiv.callback(
            exon1a="lpa_dipcn.exon1A.dipCN.txt",
            exon1b="lpa_dipcn.exon1B.dipCN.txt",
            output="lpa_kiv2_cn.txt",
            format="txt",
        )

    console.rule(f"[bold green]Pipeline complete for steps {start_step}–{end_step}[/bold green]")


if __name__ == "__main__":
    main()
