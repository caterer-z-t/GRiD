"""
Additional normalize_mosdepth tests.

Internal helper functions (process_one_individual,
compute_population_mean_depths, find_bed_gz_for_individual) were consolidated
into the top-level normalize_mosdepth() orchestrator in the current
implementation. write_read_results from count_reads remains available.

Integration-level tests are in test_pipeline.py.
"""

import threading
import pytest
from pathlib import Path

from grid.utils.count_reads import write_read_results

# --- write_read_results (count_reads) ---


def test_write_read_results(tmp_path):
    f = tmp_path / "out.tsv"
    f.write_text("")
    lock = threading.Lock()
    write_read_results(f, "S1", 123, lock)
    assert "S1\t123\n" in f.read_text()


def test_write_read_results_error(tmp_path):
    f = tmp_path / "out.tsv"
    f.write_text("")
    lock = threading.Lock()
    write_read_results(f, "S1", "Error", lock)
    assert "S1\tError\n" in f.read_text()


def test_normalize_extra_placeholder():
    """Placeholder — unit tests for process_one_individual and
    compute_population_mean_depths removed when the implementation was
    consolidated into normalize_mosdepth().
    See test_pipeline.py for integration coverage."""
    pass
