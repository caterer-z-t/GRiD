"""
Tests for compute_diploid_cn_for_exon, write_dipcn_output, and the main
compute_diploid_genotypes pipeline function (with file-based inputs).
"""
import gzip
import pytest
from pathlib import Path

from grid.utils.compute_dipcn_dir.compute_diploid_cn import compute_diploid_cn_for_exon
from grid.utils.compute_dipcn_dir.write_dipcn_output import write_dipcn_output
from grid.utils.compute_dipcn import load_neighbors


# --- compute_diploid_cn_for_exon ---

def make_counts():
    return {
        "S1": {"1B_KIV3": 100, "1B_KIV2": 50, "1B_tied": 20, "1A": 80},
        "S2": {"1B_KIV3": 90,  "1B_KIV2": 45, "1B_tied": 15, "1A": 75},
        "S3": {"1B_KIV3": 110, "1B_KIV2": 55, "1B_tied": 25, "1A": 85},
    }

def make_neighbors():
    # (scale, [(nbr_id, nbr_scale, distance), ...])
    return {
        "S1": (1.0, [("S2", 1.0, 0.1), ("S3", 1.0, 0.2)]),
        "S2": (1.0, [("S1", 1.0, 0.1), ("S3", 1.0, 0.15)]),
        "S3": (1.0, [("S1", 1.0, 0.2), ("S2", 1.0, 0.15)]),
    }

def test_compute_diploid_cn_basic():
    counts = make_counts()
    neighbors = make_neighbors()
    result = compute_diploid_cn_for_exon(counts, neighbors, "1A")
    assert set(result.keys()) == {"S1", "S2", "S3"}
    # With equal scales and similar counts, dipCN should be near 1.0
    for v in result.values():
        assert 0.5 < v < 2.0

def test_compute_diploid_cn_skips_missing_sample():
    counts = {"S1": {"1B_KIV3": 100, "1B_KIV2": 50, "1B_tied": 20, "1A": 80}}
    neighbors = {
        "S1": (1.0, [("S2", 1.0, 0.1)]),  # S2 not in counts
        "S2": (1.0, [("S1", 1.0, 0.1)]),  # S2 not in counts either
    }
    result = compute_diploid_cn_for_exon(counts, neighbors, "1A")
    # S1 has no valid neighbors → should not appear
    assert "S1" not in result

def test_compute_diploid_cn_zero_count_skipped():
    counts = {
        "S1": {"1B_KIV3": 0, "1B_KIV2": 0, "1B_tied": 0, "1A": 0},
        "S2": {"1B_KIV3": 50, "1B_KIV2": 20, "1B_tied": 5, "1A": 40},
    }
    neighbors = {"S1": (1.0, [("S2", 1.0, 0.1)]), "S2": (1.0, [("S1", 1.0, 0.1)])}
    result = compute_diploid_cn_for_exon(counts, neighbors, "1A")
    assert "S1" not in result  # zero count → skipped

def test_compute_diploid_cn_n_neighbors_limit():
    counts = make_counts()
    neighbors = make_neighbors()
    # n_neighbors=1 → only use first neighbor
    result_1 = compute_diploid_cn_for_exon(counts, neighbors, "1A", n_neighbors=1)
    result_all = compute_diploid_cn_for_exon(counts, neighbors, "1A", n_neighbors=200)
    # Both should produce results; values may differ
    assert set(result_1.keys()) == set(result_all.keys())

def test_compute_diploid_cn_exon_types():
    counts = make_counts()
    neighbors = make_neighbors()
    for exon in ("1B_KIV3", "1B_notKIV3", "1B", "1A"):
        result = compute_diploid_cn_for_exon(counts, neighbors, exon)
        assert len(result) > 0


# --- write_dipcn_output ---

def test_write_dipcn_output(tmp_path):
    results = {"S2": 2.5, "S1": 1.8, "S3": 3.1}
    out = tmp_path / "dipcn.tsv"
    write_dipcn_output(results, str(out))
    lines = out.read_text().strip().split("\n")
    assert lines[0] == "ID\tdipCN"
    # Should be sorted by sample ID
    ids = [l.split("\t")[0] for l in lines[1:]]
    assert ids == sorted(results.keys())

def test_write_dipcn_output_creates_dirs(tmp_path):
    out = tmp_path / "deep" / "dir" / "dipcn.tsv"
    write_dipcn_output({"S1": 2.0}, str(out))
    assert out.exists()

def test_write_dipcn_output_format(tmp_path):
    out = tmp_path / "dipcn.tsv"
    write_dipcn_output({"S1": 2.123456789}, str(out))
    line = out.read_text().strip().split("\n")[1]
    _, val = line.split("\t")
    # Should be written to 6 decimal places
    assert val == "2.123457"


# --- compute_diploid_genotypes (main function) via temp files ---

def test_compute_diploid_genotypes_end_to_end(tmp_path):
    """Integration test: build minimal read-counts + neighbors files and run compute_diploid_genotypes."""
    from contextlib import contextmanager
    from unittest.mock import patch, MagicMock
    from grid.utils.compute_dipcn import compute_diploid_genotypes

    # Write read counts file
    counts_file = tmp_path / "counts.tsv"
    counts_file.write_text("Sample\tchr6:1000-2000\nS1\t100\nS2\t90\nS3\t110\n")

    # Write neighbors file
    nbrs_file = tmp_path / "nbrs.zMax2.0.tsv.gz"
    with gzip.open(nbrs_file, "wt") as f:
        f.write("S1\t1.0\tS2\t1.0\t0.1\tS3\t1.0\t0.2\n")
        f.write("S2\t1.0\tS1\t1.0\t0.1\tS3\t1.0\t0.15\n")
        f.write("S3\t1.0\tS1\t1.0\t0.2\tS2\t1.0\t0.15\n")

    config = {
        "output_dir": str(tmp_path),
        "output_file_type": "tsv",
        "compute_diploid_genotypes": {"output_file_prefix": "dipcn", "n_nbr": 10},
        "count_reads": {"output_file_prefix": "counts"},
        "mosdepth": {"neighbors": {"zmax": 2.0, "output_file_prefix": "nbrs"}},
    }

    # progress_bar uses custom Rich theme styles only available on the CLI console;
    # stub it out so the test works without a themed console.
    @contextmanager
    def _noop_progress(*args, **kwargs):
        mock_progress = MagicMock()
        mock_task = MagicMock()
        yield mock_progress, mock_task

    with patch("grid.utils.compute_dipcn.progress_bar", _noop_progress):
        compute_diploid_genotypes(config, console=None)

    out = tmp_path / "dipcn.tsv"
    assert out.exists()
    lines = out.read_text().strip().split("\n")
    # Header + 3 data rows
    assert len(lines) == 4
    assert lines[0] == "Sample\tNorm_Reads"
