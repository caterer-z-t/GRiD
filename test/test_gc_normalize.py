"""Tests for gc_normalize.py and compute_gc_content.py."""
import gzip
import math
import tempfile
from pathlib import Path

import numpy as np
import pytest

from grid.utils.gc_normalize import (
    _lowess_fit,
    _lowess_predict,
    _mad,
    _read_depths_from_mosdepth,
)
from grid.utils.compute_gc_content import load_gc_content_file, write_gc_content_file
from grid.utils.select_norm_regions import _filter_regions, _select_gc_diverse


# ---------------------------------------------------------------------------
# _mad
# ---------------------------------------------------------------------------

def test_mad_simple():
    arr = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    assert abs(_mad(arr) - 1.0) < 1e-9


def test_mad_constant():
    arr = np.array([5.0] * 10)
    assert _mad(arr) == 0.0


def test_mad_empty():
    assert math.isnan(_mad(np.array([])))


# ---------------------------------------------------------------------------
# _lowess_fit / _lowess_predict
# ---------------------------------------------------------------------------

@pytest.fixture
def linear_gc_depth():
    """Synthetic linear relationship: depth = 30 + 20 * gc + noise."""
    rng = np.random.default_rng(0)
    gcs = np.linspace(0.2, 0.85, 300)
    depths = 30.0 + 20.0 * gcs + rng.normal(0, 0.3, 300)
    return gcs, depths


def test_lowess_fit_returns_grid(linear_gc_depth):
    gcs, depths = linear_gc_depth
    grid_x, grid_y = _lowess_fit(gcs, depths, frac=0.3, n_grid=50)
    assert len(grid_x) == 50
    assert len(grid_y) == 50
    assert np.all(np.diff(grid_x) > 0), "grid_x should be strictly increasing"


def test_lowess_predicts_linear_trend(linear_gc_depth):
    gcs, depths = linear_gc_depth
    grid_x, grid_y = _lowess_fit(gcs, depths, frac=0.3, n_grid=50)
    # Predict at gc=0.65: expected 30 + 20*0.65 = 43.0
    pred = float(np.interp(0.65, grid_x, grid_y))
    assert abs(pred - 43.0) < 2.0, f"Prediction: {pred:.2f}"


def test_lowess_predict_wrapper(linear_gc_depth):
    gcs, depths = linear_gc_depth
    pred = _lowess_predict(gcs, depths, 0.65, frac=0.3)
    assert pred is not None
    assert abs(pred - 43.0) < 2.0


def test_lowess_predict_constant_gc():
    # Edge case: all GC values are the same
    gcs = np.full(50, 0.5)
    depths = np.full(50, 30.0)
    grid_x, grid_y = _lowess_fit(gcs, depths, frac=0.3)
    # Should not raise and grid should have correct shape
    assert len(grid_x) == len(grid_y)


# ---------------------------------------------------------------------------
# _read_depths_from_mosdepth
# ---------------------------------------------------------------------------

def _make_mosdepth_gz(tmp_path: Path, rows: list[tuple]) -> Path:
    """Create a fake mosdepth regions.bed.gz for testing."""
    bed_gz = tmp_path / "sample.regions.bed.gz"
    with gzip.open(bed_gz, "wt") as f:
        for chrom, start, end, depth in rows:
            f.write(f"{chrom}\t{start}\t{end}\t{depth}\n")
    return bed_gz


def test_read_depths_norm_and_kiv2(tmp_path):
    gc_map = {
        ("chr1", 0, 1000): 0.42,
        ("chr1", 1000, 2000): 0.55,
    }
    rows = [
        ("chr1", 0, 1000, 30.0),
        ("chr1", 1000, 2000, 35.0),
        ("chr6", 160611568, 160646868, 28.0),  # KIV-2 target
    ]
    bed_gz = _make_mosdepth_gz(tmp_path, rows)
    norm_depths, kiv2_depth = _read_depths_from_mosdepth(
        bed_gz, gc_map,
        target_chrom="chr6",
        target_start=160611568,
        target_end=160646868,
    )
    assert len(norm_depths) == 2
    assert kiv2_depth is not None
    assert abs(kiv2_depth - 28.0) < 1e-6


def test_read_depths_no_target_overlap(tmp_path):
    """If no bins overlap the target region, kiv2_depth should be None."""
    gc_map = {("chr1", 0, 1000): 0.42}
    rows = [("chr1", 0, 1000, 30.0)]
    bed_gz = _make_mosdepth_gz(tmp_path, rows)
    _, kiv2_depth = _read_depths_from_mosdepth(
        bed_gz, gc_map,
        target_chrom="chr6",
        target_start=160611568,
        target_end=160646868,
    )
    assert kiv2_depth is None


# ---------------------------------------------------------------------------
# GC content file I/O
# ---------------------------------------------------------------------------

def test_write_load_gc_content_file(tmp_path):
    gc_map = {
        ("chr1", 0, 1000): 0.45,
        ("chr2", 5000, 6000): 0.62,
    }
    out = tmp_path / "gc.tsv"
    write_gc_content_file(gc_map, out)
    loaded = load_gc_content_file(out)
    assert len(loaded) == 2
    assert abs(loaded[("chr1", 0, 1000)] - 0.45) < 1e-5
    assert abs(loaded[("chr2", 5000, 6000)] - 0.62) < 1e-5


# ---------------------------------------------------------------------------
# _filter_regions
# ---------------------------------------------------------------------------

def test_filter_regions_basic():
    sums = {
        ("chr1", 0, 1000): 300.0,        # mean=30, CV~0.01 → passes
        ("chr1", 1000, 2000): 3000.0,    # mean=300, > max_depth=100 → fails
        ("chr1", 2000, 3000): 25.0,      # mean=2.5, < min_depth=20 → fails
    }
    sq_s = {
        ("chr1", 0, 1000): 9001.0,
        ("chr1", 1000, 2000): 90001.0,
        ("chr1", 2000, 3000): 63.0,
    }
    cnts = {k: 10 for k in sums}
    passing = _filter_regions(sums, sq_s, cnts, min_depth=20, max_depth=100, max_cv=0.15)
    assert len(passing) == 1
    assert passing[0][0] == ("chr1", 0, 1000)


def test_filter_regions_high_cv_excluded():
    # mean=30, std=10 → CV=0.33 > 0.15
    n = 20
    sums = {("chr1", 0, 1000): 30.0 * n}
    sq_s = {("chr1", 0, 1000): (30.0 ** 2 + 10.0 ** 2) * n}
    cnts = {("chr1", 0, 1000): n}
    passing = _filter_regions(sums, sq_s, cnts, min_depth=20, max_depth=100, max_cv=0.15)
    assert len(passing) == 0


# ---------------------------------------------------------------------------
# _select_gc_diverse
# ---------------------------------------------------------------------------

def test_select_gc_diverse_one_per_bin():
    regions_cv = [
        (("chr1", 0, 1000), 0.05),
        (("chr1", 1000, 2000), 0.08),
        (("chr1", 2000, 3000), 0.03),
        (("chr1", 3000, 4000), 0.10),
    ]
    gc_map = {
        ("chr1", 0, 1000): 0.42,
        ("chr1", 1000, 2000): 0.43,   # same GC bin as above
        ("chr1", 2000, 3000): 0.67,
        ("chr1", 3000, 4000): 0.68,   # same GC bin as above
    }
    selected = _select_gc_diverse(regions_cv, gc_map, gc_bin_size=0.05, n_per_gc_bin=1)
    assert len(selected) == 2
    assert ("chr1", 0, 1000) in selected    # CV=0.05, lowest in its bin
    assert ("chr1", 2000, 3000) in selected  # CV=0.03, lowest in its bin


def test_select_gc_diverse_nan_gc_excluded():
    regions_cv = [(("chr1", 0, 1000), 0.05)]
    gc_map = {("chr1", 0, 1000): math.nan}
    selected = _select_gc_diverse(regions_cv, gc_map, gc_bin_size=0.05, n_per_gc_bin=10)
    assert len(selected) == 0
