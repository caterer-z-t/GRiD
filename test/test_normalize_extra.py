"""
Additional normalize_mosdepth tests covering process_one_individual,
compute_population_mean_depths, and write_read_results from count_reads.
"""
import gzip
import pytest
from pathlib import Path

from grid.utils.normalize_mosdepth import (
    process_one_individual,
    compute_population_mean_depths,
    find_bed_gz_for_individual,
)
from grid.utils.count_reads import write_read_results


# ── helpers ────────────────────────────────────────────────────────────────

def write_bed_gz(path, rows):
    """Write a minimal mosdepth regions.bed.gz file."""
    with gzip.open(path, "wt") as f:
        for chrom, start, end, depth in rows:
            f.write(f"{chrom}\t{start}\t{end}\t{depth}\n")


# ── process_one_individual ─────────────────────────────────────────────────

def test_process_one_individual_basic(tmp_path):
    bed = tmp_path / "S1.regions.bed.gz"
    write_bed_gz(bed, [
        ("chr6", 1000, 2000, 30.0),
        ("chr6", 2000, 3000, 40.0),
    ])
    valid = {(1000, 2000), (2000, 3000)}
    ind_id, results = process_one_individual("S1", str(tmp_path), "chr6", 1000, 3000, valid, {})
    assert ind_id == "S1"
    assert len(results) == 2

def test_process_one_individual_filters_out_of_range(tmp_path):
    bed = tmp_path / "S1.regions.bed.gz"
    write_bed_gz(bed, [
        ("chr6", 5000, 6000, 30.0),  # outside 1000–3000
        ("chr6", 1000, 2000, 25.0),
    ])
    valid = {(1000, 2000)}
    _, results = process_one_individual("S1", str(tmp_path), "chr6", 1000, 3000, valid, {})
    assert len(results) == 1

def test_process_one_individual_filters_not_in_valid_regions(tmp_path):
    bed = tmp_path / "S1.regions.bed.gz"
    write_bed_gz(bed, [("chr6", 1000, 2000, 30.0)])
    valid = {(2000, 3000)}  # (1000,2000) not valid
    _, results = process_one_individual("S1", str(tmp_path), "chr6", 1000, 3000, valid, {})
    assert results == []

def test_process_one_individual_respects_repeat_mask(tmp_path):
    bed = tmp_path / "S1.regions.bed.gz"
    write_bed_gz(bed, [("chr6", 1000, 2000, 30.0)])
    valid = {(1000, 2000)}
    # Exclude kb bin 1 (covers 1000–1999)
    excluded = {"chr6": {1}}
    _, results = process_one_individual("S1", str(tmp_path), "chr6", 1000, 3000, valid, excluded)
    assert results == []

def test_process_one_individual_missing_file(tmp_path):
    valid = {(1000, 2000)}
    ind_id, results = process_one_individual("NOSUCH", str(tmp_path), "chr6", 1000, 3000, valid, {})
    assert ind_id == "NOSUCH"
    assert results == []

def test_process_one_individual_wrong_chrom(tmp_path):
    bed = tmp_path / "S1.regions.bed.gz"
    write_bed_gz(bed, [("chr1", 1000, 2000, 30.0)])
    valid = {(1000, 2000)}
    _, results = process_one_individual("S1", str(tmp_path), "chr6", 1000, 3000, valid, {})
    assert results == []

def test_process_one_individual_zero_depth_filtered(tmp_path):
    bed = tmp_path / "S1.regions.bed.gz"
    write_bed_gz(bed, [("chr6", 1000, 2000, 0.0)])
    valid = {(1000, 2000)}
    _, results = process_one_individual("S1", str(tmp_path), "chr6", 1000, 3000, valid, {})
    assert results == []


# ── compute_population_mean_depths ────────────────────────────────────────

def test_compute_population_mean_depths_basic(tmp_path):
    for sid in ["S1", "S2"]:
        bed = tmp_path / f"{sid}.regions.bed.gz"
        write_bed_gz(bed, [
            ("chr6", 1000, 2000, 30.0),
            ("chr6", 2000, 3000, 40.0),
        ])
    individuals = {
        "S1": tmp_path / "S1.regions.bed.gz",
        "S2": tmp_path / "S2.regions.bed.gz",
    }
    means = compute_population_mean_depths(individuals, str(tmp_path), "chr6", 1000, 3000, {}, threads=1)
    assert (1000, 2000) in means
    assert means[(1000, 2000)] == pytest.approx(30.0)
    assert means[(2000, 3000)] == pytest.approx(40.0)

def test_compute_population_mean_depths_excludes_masked_regions(tmp_path):
    bed = tmp_path / "S1.regions.bed.gz"
    write_bed_gz(bed, [("chr6", 1000, 2000, 30.0)])
    individuals = {"S1": bed}
    excluded = {"chr6": {1}}  # kb bin 1 → covers 1000–1999
    means = compute_population_mean_depths(individuals, str(tmp_path), "chr6", 1000, 3000, excluded, threads=1)
    assert (1000, 2000) not in means

def test_compute_population_mean_depths_missing_file(tmp_path):
    individuals = {"NOSUCH": tmp_path / "NOSUCH.regions.bed.gz"}
    means = compute_population_mean_depths(individuals, str(tmp_path), "chr6", 1000, 3000, {}, threads=1)
    assert means == {}


# ── write_read_results (count_reads) ──────────────────────────────────────

def test_write_read_results(tmp_path):
    import threading
    f = tmp_path / "out.tsv"
    f.write_text("")
    lock = threading.Lock()
    write_read_results(f, "S1", 123, lock)
    assert "S1\t123\n" in f.read_text()

def test_write_read_results_error(tmp_path):
    import threading
    f = tmp_path / "out.tsv"
    f.write_text("")
    lock = threading.Lock()
    write_read_results(f, "S1", "Error", lock)
    assert "S1\tError\n" in f.read_text()
