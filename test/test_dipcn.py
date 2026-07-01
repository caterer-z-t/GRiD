import gzip
import pytest
from pathlib import Path

from grid.utils.compute_dipcn import load_neighbors
from grid.utils.compute_dipcn_dir.get_exon_count import get_exon_count
from grid.utils.compute_dipcn_dir.normalize_sample_id import normalize_sample_id
from grid.utils.compute_dipcn_dir.validate_sample_overlap import validate_sample_overlap
from grid.utils.compute_dipcn_dir.load_count_results import load_count_results
from grid.utils.compute_dipcn_dir.load_neighbor_results import load_neighbor_results


# --- normalize_sample_id ---

def test_normalize_sample_id_plain():
    assert normalize_sample_id("NWD123") == "NWD123"

def test_normalize_sample_id_cram():
    assert normalize_sample_id("NWD123.cram") == "NWD123"

def test_normalize_sample_id_bam():
    assert normalize_sample_id("NWD123.bam") == "NWD123"

def test_normalize_sample_id_b38_suffix():
    assert normalize_sample_id("NWD123.b38.irc.v1_subset") == "NWD123"

def test_normalize_sample_id_strips_whitespace():
    assert normalize_sample_id("  NWD123  ") == "NWD123"


# --- get_exon_count ---

def test_get_exon_count_1B_KIV3():
    counts = {"1B_KIV3": 10, "1B_KIV2": 5, "1B_tied": 3, "1A": 7}
    assert get_exon_count(counts, "1B_KIV3") == 10

def test_get_exon_count_1B_notKIV3():
    counts = {"1B_KIV3": 10, "1B_KIV2": 5, "1B_tied": 3, "1A": 7}
    assert get_exon_count(counts, "1B_notKIV3") == 8  # 5+3

def test_get_exon_count_1B():
    counts = {"1B_KIV3": 10, "1B_KIV2": 5, "1B_tied": 3, "1A": 7}
    assert get_exon_count(counts, "1B") == 18  # 10+5+3

def test_get_exon_count_1A():
    counts = {"1B_KIV3": 10, "1B_KIV2": 5, "1B_tied": 3, "1A": 7}
    assert get_exon_count(counts, "1A") == 7

def test_get_exon_count_unknown():
    with pytest.raises(ValueError):
        get_exon_count({}, "BADTYPE")

def test_get_exon_count_missing_keys():
    # missing keys default to 0
    assert get_exon_count({}, "1B") == 0


# --- validate_sample_overlap ---

def test_validate_sample_overlap_full():
    counts = {"S1": {}, "S2": {}}
    neighbors = {"S1": None, "S2": None}
    n, overlap = validate_sample_overlap(counts, neighbors)
    assert n == 2
    assert overlap == {"S1", "S2"}

def test_validate_sample_overlap_partial():
    counts = {"S1": {}, "S2": {}}
    neighbors = {"S2": None, "S3": None}
    n, overlap = validate_sample_overlap(counts, neighbors)
    assert n == 1
    assert overlap == {"S2"}

def test_validate_sample_overlap_none():
    n, overlap = validate_sample_overlap({"A": {}}, {"B": None})
    assert n == 0
    assert overlap == set()


# --- load_count_results ---

def test_load_count_results(tmp_path):
    f = tmp_path / "counts.tsv"
    f.write_text("NWD123\t10\t5\t3\t7\nNWD456\t2\t1\t0\t4\n")
    result = load_count_results(f)
    assert "NWD123" in result
    assert result["NWD123"]["1B_KIV3"] == 10
    assert result["NWD456"]["1A"] == 4

def test_load_count_results_skips_bad_lines(tmp_path):
    f = tmp_path / "counts.tsv"
    f.write_text("NWD123\t10\t5\t3\t7\nbadline\nNWD456\t2\t1\t0\t4\n")
    result = load_count_results(f)
    assert len(result) == 2

def test_load_count_results_normalizes_ids(tmp_path):
    f = tmp_path / "counts.tsv"
    f.write_text("NWD123.cram\t10\t5\t3\t7\n")
    result = load_count_results(f)
    assert "NWD123" in result


# --- load_neighbors (compute_dipcn) ---

def test_load_neighbors(tmp_path):
    nbr_file = tmp_path / "nbrs.tsv.gz"
    with gzip.open(nbr_file, "wt") as f:
        # format: id  scale  nbr1_id  nbr1_scale  dist  nbr2_id ...
        f.write("S1\t30.5\tS2\t28.0\t0.05\tS3\t31.0\t0.08\n")
        f.write("S2\t28.0\tS1\t30.5\t0.05\n")

    neighbors, scales = load_neighbors(nbr_file)
    assert "S1" in neighbors
    assert scales["S1"] == pytest.approx(30.5)
    assert len(neighbors["S1"]) == 2
    # Each entry is (neighbor_id, neighbor_scale)
    nbr_ids = [n for n, _ in neighbors["S1"]]
    assert "S2" in nbr_ids

def test_load_neighbors_empty(tmp_path):
    nbr_file = tmp_path / "nbrs.tsv.gz"
    with gzip.open(nbr_file, "wt") as f:
        f.write("")
    neighbors, scales = load_neighbors(nbr_file)
    assert neighbors == {}
    assert scales == {}


# --- load_neighbor_results (compute_dipcn_dir) ---

def test_load_neighbor_results(tmp_path):
    f = tmp_path / "nbrs.tsv.gz"
    with gzip.open(f, "wt") as fh:
        fh.write("NWD123\t1.0\tNWD456\t0.98\t0.05\tNWD789\t1.02\t0.07\n")
    result = load_neighbor_results(f)
    assert "NWD123" in result
    scale, nbr_list = result["NWD123"]
    assert scale == pytest.approx(1.0)
    assert len(nbr_list) == 2
    assert nbr_list[0][0] == "NWD456"
