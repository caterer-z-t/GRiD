"""
mosdepth tests.

Internal helper functions (compute_region_coverage, build_mosdepth_command,
check_mosdepth_available, write_coverage_result) were consolidated into
compute_mosdepth / run_mosdepth_single_cram in the current implementation.
The tests below cover the remaining public symbol: remove_intermediate_files.
Integration-level tests for compute_mosdepth are in test_pipeline.py.
"""

import pytest
from pathlib import Path

from grid.utils.mosdepth import remove_intermediate_files

# --- remove_intermediate_files ---


def test_remove_intermediate_files(tmp_path):
    (tmp_path / "sample.mosdepth.global.dist.txt").touch()
    (tmp_path / "sample.regions.bed.gz.csi").touch()
    (tmp_path / "keep_me.tsv").touch()
    remove_intermediate_files(tmp_path)
    assert not (tmp_path / "sample.mosdepth.global.dist.txt").exists()
    assert not (tmp_path / "sample.regions.bed.gz.csi").exists()
    assert (tmp_path / "keep_me.tsv").exists()


def test_remove_intermediate_files_removes_bed_gz(tmp_path):
    (tmp_path / "sample.regions.bed.gz").touch()
    remove_intermediate_files(tmp_path)
    assert not (tmp_path / "sample.regions.bed.gz").exists()


def test_remove_intermediate_files_empty_dir(tmp_path):
    remove_intermediate_files(tmp_path)  # should not raise
