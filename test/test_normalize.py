import gzip
import numpy as np
import pytest
from pathlib import Path

from grid.utils.normalize_mosdepth import (
    norm_chrom,
    load_repeat_mask,
    normalize_matrix,
    select_high_variance_regions,
    build_matrix_from_regions,
    filter_empty_samples,
    map_mosdepth_files_to_samples,
    find_bed_gz_for_individual,
    write_normalized_output,
)
from grid.utils.find_neighbors import read_normalized_data


# --- norm_chrom ---

def test_norm_chrom_already_prefixed():
    assert norm_chrom("chr6") == "chr6"

def test_norm_chrom_adds_prefix():
    assert norm_chrom("6") == "chr6"

def test_norm_chrom_x():
    assert norm_chrom("X") == "chrX"


# --- load_repeat_mask ---

def test_load_repeat_mask(tmp_path):
    bed = tmp_path / "mask.bed"
    bed.write_text("chr6\t1000\t3000\n")
    excluded = load_repeat_mask(str(bed))
    assert "chr6" in excluded
    # bins 1–3 (1000//1000=1, 3000//1000=3)
    assert 1 in excluded["chr6"]
    assert 3 in excluded["chr6"]

def test_load_repeat_mask_skips_comments(tmp_path):
    bed = tmp_path / "mask.bed"
    bed.write_text("# comment\nchr1\t0\t1000\n")
    excluded = load_repeat_mask(str(bed))
    assert "chr1" in excluded

def test_load_repeat_mask_normalizes_chrom(tmp_path):
    bed = tmp_path / "mask.bed"
    bed.write_text("6\t0\t1000\n")
    excluded = load_repeat_mask(str(bed))
    assert "chr6" in excluded

def test_load_repeat_mask_skips_short_lines(tmp_path):
    bed = tmp_path / "mask.bed"
    bed.write_text("chr1\t0\n")  # only 2 fields
    excluded = load_repeat_mask(str(bed))
    assert excluded == {}


# --- build_matrix_from_regions ---

def test_build_matrix_basic():
    regions = {
        "S1": [(0, 1000, 30.0), (1000, 2000, 40.0)],
        "S2": [(0, 1000, 20.0), (1000, 2000, 25.0)],
    }
    order, mat = build_matrix_from_regions(regions)
    assert mat.shape == (2, 2)
    assert set(order) == {"S1", "S2"}

def test_build_matrix_missing_region():
    regions = {
        "S1": [(0, 1000, 30.0)],
        "S2": [(0, 1000, 20.0), (1000, 2000, 25.0)],
    }
    order, mat = build_matrix_from_regions(regions)
    assert mat.shape == (2, 2)
    # S1 is missing region (1000,2000) → NaN
    i = order.index("S1")
    assert np.isnan(mat[i, 1])


# --- normalize_matrix ---

def test_normalize_matrix_shape():
    mat = np.array([[30.0, 40.0, 35.0],
                    [20.0, 25.0, 22.0],
                    [35.0, 45.0, 40.0]])
    norm, ratios, col_means, col_vars = normalize_matrix(mat)
    assert norm.shape == mat.shape

def test_normalize_matrix_returns_variance_ratios():
    mat = np.array([[30.0, 40.0],
                    [20.0, 60.0],
                    [40.0, 20.0]])
    _, ratios, _, _ = normalize_matrix(mat)
    assert isinstance(ratios, dict)
    assert len(ratios) == 2


# --- select_high_variance_regions ---

def test_select_high_variance_regions():
    ratios = {0: 1.0, 1: 5.0, 2: 10.0, 3: 2.0}
    # top_frac=0.5 → keep top 50% by variance
    selected = select_high_variance_regions(ratios, top_frac=0.5)
    # threshold is the value at index int(0.5*4)=2 of sorted [1,2,5,10] → 5
    # keep indices where ratio > 5 → only index 2 (ratio=10)
    assert 2 in selected

def test_select_high_variance_regions_empty():
    assert select_high_variance_regions({}) == []


# --- filter_empty_samples ---

def test_filter_empty_samples_removes_empty():
    data = {"S1": [(0, 1000, 30.0)], "S2": [], "S3": [(0, 1000, 20.0)]}
    filtered = filter_empty_samples(data)
    assert "S2" not in filtered
    assert "S1" in filtered
    assert "S3" in filtered

def test_filter_empty_samples_all_empty():
    filtered = filter_empty_samples({"S1": [], "S2": []})
    assert filtered == {}


# --- map_mosdepth_files_to_samples ---

def test_map_mosdepth_files_to_samples(tmp_path):
    (tmp_path / "S1.regions.bed.gz").touch()
    (tmp_path / "S2.regions.bed.gz").touch()
    (tmp_path / "other.txt").touch()
    result = map_mosdepth_files_to_samples(str(tmp_path), ["S1", "S2", "S3"])
    assert "S1" in result
    assert "S2" in result
    assert "S3" not in result  # no file for S3


# --- find_bed_gz_for_individual ---

def test_find_bed_gz_for_individual_found(tmp_path):
    f = tmp_path / "S1.regions.bed.gz"
    f.touch()
    result = find_bed_gz_for_individual("S1", str(tmp_path))
    assert result.exists()

def test_find_bed_gz_for_individual_not_found(tmp_path):
    result = find_bed_gz_for_individual("NOSUCH", str(tmp_path))
    assert not result.exists()


# --- write_normalized_output / read_normalized_data round-trip ---

def test_normalize_write_read_roundtrip(tmp_path):
    mat = np.array([[1.0, 2.0, 3.0],
                    [4.0, 5.0, 6.0]])
    individuals = ["S1", "S2"]
    selected = [0, 2]
    col_means = np.array([2.5, 3.5, 4.5])
    col_vars = np.array([4.5, 4.5, 4.5])
    raw_means = np.array([2.0, 5.0])

    out = tmp_path / "norm.tsv.gz"
    write_normalized_output(mat, individuals, selected, out, col_means, col_vars, raw_means)

    assert out.exists()
    indivs, sigma2ratios, data_matrix, scales = read_normalized_data(out)

    assert indivs == ["S1", "S2"]
    assert data_matrix.shape == (2, 2)  # 2 samples × 2 selected regions
    assert scales["S1"] == pytest.approx(2.0, abs=0.01)
    assert scales["S2"] == pytest.approx(5.0, abs=0.01)
