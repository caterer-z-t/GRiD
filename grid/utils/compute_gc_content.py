"""
Compute GC content for genomic regions from a reference FASTA.
"""
from __future__ import annotations

import gzip
import math
from pathlib import Path
from typing import Iterable

import pysam


def compute_gc_for_region(fasta: pysam.FastaFile, chrom: str, start: int, end: int) -> float:
    """
    Return GC fraction [0.0, 1.0] for reference region [start, end).
    Returns NaN if the region is empty or contains only Ns.
    """
    try:
        seq = fasta.fetch(chrom, start, end).upper()
    except (ValueError, KeyError):
        return math.nan
    if not seq:
        return math.nan
    atgc = sum(seq.count(b) for b in "ATGC")
    if atgc == 0:
        return math.nan
    return (seq.count("G") + seq.count("C")) / atgc


def compute_gc_for_regions(
    fasta_path: str | Path,
    regions: Iterable[tuple[str, int, int]],
) -> dict[tuple[str, int, int], float]:
    """
    Compute GC content for an iterable of (chrom, start, end) tuples.

    Returns:
        dict mapping (chrom, start, end) → gc_fraction
    """
    fasta_path = str(fasta_path)
    gc_map: dict[tuple[str, int, int], float] = {}
    with pysam.FastaFile(fasta_path) as fasta:
        for chrom, start, end in regions:
            gc_map[(chrom, start, end)] = compute_gc_for_region(fasta, chrom, start, end)
    return gc_map


def load_gc_content_file(
    gc_file: str | Path,
) -> dict[tuple[str, int, int], float]:
    """
    Load a pre-computed GC content file produced by select_norm_regions.

    Expected format (tab-separated, no header):
        chrom  start  end  gc_fraction

    Returns:
        dict mapping (chrom, start, end) → gc_fraction
    """
    gc_file = Path(gc_file)
    gc_map: dict[tuple[str, int, int], float] = {}
    opener = gzip.open if str(gc_file).endswith(".gz") else open
    with opener(gc_file, "rt") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 4:
                continue
            try:
                chrom = parts[0]
                start = int(parts[1])
                end = int(parts[2])
                gc = float(parts[3])
                gc_map[(chrom, start, end)] = gc
            except ValueError:
                continue
    return gc_map


def write_gc_content_file(
    gc_map: dict[tuple[str, int, int], float],
    output_file: str | Path,
) -> None:
    """Write GC content map to a tab-separated file."""
    output_file = Path(output_file)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    with open(output_file, "w") as out:
        for (chrom, start, end), gc in sorted(gc_map.items()):
            out.write(f"{chrom}\t{start}\t{end}\t{gc:.6f}\n")
