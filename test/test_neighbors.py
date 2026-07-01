import gzip
import numpy as np
import pytest
from pathlib import Path

from grid.utils.find_neighbors import (
    filter_regions_by_variance,
    find_neighbors_sklearn,
    save_neighbors,
    read_normalized_data,
)
from grid.utils.normalize_mosdepth import write_normalized_output


# --- filter_regions_by_variance ---

def test_filter_regions_basic():
    ratios = np.array([1.0, 5.0, 1000.5, 200.0])
    valid, R_use = filter_regions_by_variance(ratios, frac_r=1.0, sigma2_max=1000.0)
    # ratio 1000.5 exceeds sigma2_max → excluded
    assert 2 not in valid
    assert R_use == len(valid)

def test_filter_regions_all_valid():
    ratios = np.array([1.0, 2.0, 3.0])
    valid, R_use = filter_regions_by_variance(ratios, frac_r=1.0, sigma2_max=1000.0)
    assert R_use == 3

def test_filter_regions_empty():
    ratios = np.array([np.nan, np.nan])
    valid, R_use = filter_regions_by_variance(ratios)
    # no finite values → keep all
    assert R_use == 2

def test_filter_regions_frac_r():
    # frac_r=0.5: lower_idx = int(4 * 0.5) = 2; sorted=[1,2,3,4] → sigma2_min=3
    ratios = np.array([1.0, 2.0, 3.0, 4.0])
    valid, R_use = filter_regions_by_variance(ratios, frac_r=0.5, sigma2_max=1000.0)
    assert R_use < 4  # some filtered by lower bound


# --- find_neighbors_sklearn ---

def test_find_neighbors_sklearn_basic():
    data = np.array([[1.0, 0.0],
                     [1.1, 0.1],
                     [5.0, 5.0]])
    individuals = ["S1", "S2", "S3"]
    result = find_neighbors_sklearn(data, individuals, n_neighbors=2)
    assert set(result.keys()) == {"S1", "S2", "S3"}
    # S1 and S2 should be each other's nearest neighbor
    s1_nbrs = [nbr for nbr, _ in result["S1"]]
    assert "S2" in s1_nbrs

def test_find_neighbors_sklearn_excludes_self():
    data = np.array([[1.0, 0.0], [2.0, 0.0], [3.0, 0.0]])
    individuals = ["A", "B", "C"]
    result = find_neighbors_sklearn(data, individuals, n_neighbors=2)
    for ind, nbrs in result.items():
        assert ind not in [n for n, _ in nbrs]

def test_find_neighbors_sklearn_fewer_than_requested():
    data = np.array([[1.0], [2.0]])
    individuals = ["A", "B"]
    result = find_neighbors_sklearn(data, individuals, n_neighbors=10)
    # only 1 neighbor available per sample
    assert len(result["A"]) == 1

def test_find_neighbors_distances_are_squared():
    data = np.array([[0.0, 0.0], [3.0, 4.0]])  # Euclidean dist = 5
    individuals = ["A", "B"]
    result = find_neighbors_sklearn(data, individuals, n_neighbors=1)
    _, sq_dist = result["A"][0]
    assert sq_dist == pytest.approx(25.0, rel=1e-5)


# --- save_neighbors / read_normalized_data round-trip via file ---

def test_save_neighbors_format(tmp_path):
    neighbors = {
        "S1": [("S2", 0.5), ("S3", 1.2)],
        "S2": [("S1", 0.5)],
    }
    scales = {"S1": 30.0, "S2": 28.0, "S3": 32.0}
    out = tmp_path / "nbrs.tsv.gz"
    save_neighbors(neighbors, scales, out, zmax=2.0, R_use=100)

    assert out.exists()
    with gzip.open(out, "rt") as f:
        lines = f.readlines()
    assert len(lines) == 2

    # First line: S1  scale  S2  scale  dist  S3  scale  dist
    fields = lines[0].strip().split("\t")
    assert fields[0] == "S1"
    assert float(fields[1]) == pytest.approx(30.0)
    # neighbor entries come in triplets starting at index 2
    assert fields[2] == "S2"

def test_save_neighbors_distance_normalization(tmp_path):
    neighbors = {"A": [("B", 200.0)]}  # sq_dist=200, R_use=100 → norm=1.0
    scales = {"A": 1.0, "B": 1.0}
    out = tmp_path / "nbrs.tsv.gz"
    save_neighbors(neighbors, scales, out, zmax=2.0, R_use=100)
    with gzip.open(out, "rt") as f:
        line = f.readline()
    fields = line.strip().split("\t")
    norm_dist = float(fields[4])
    assert norm_dist == pytest.approx(1.0, rel=1e-5)

def test_save_neighbors_r_use_zero_guard(tmp_path):
    # R_use=0 should not cause divide-by-zero
    neighbors = {"A": [("B", 10.0)]}
    scales = {"A": 1.0, "B": 1.0}
    out = tmp_path / "nbrs.tsv.gz"
    save_neighbors(neighbors, scales, out, zmax=2.0, R_use=0)
    assert out.exists()


# --- read_normalized_data ---

def test_read_normalized_data_roundtrip(tmp_path):
    mat = np.array([[0.5, -0.5], [-0.5, 0.5]])
    individuals = ["X1", "X2"]
    col_means = np.array([1.0, 2.0])
    col_vars = np.array([0.1, 0.2])
    raw_means = np.array([25.0, 30.0])
    out = tmp_path / "norm.tsv.gz"

    write_normalized_output(mat, individuals, [0, 1], out, col_means, col_vars, raw_means)
    indivs, ratios, data, scales = read_normalized_data(out)

    assert indivs == ["X1", "X2"]
    assert data.shape == (2, 2)
    assert scales["X1"] == pytest.approx(25.0, abs=0.01)
    assert ratios.shape == (2,)
