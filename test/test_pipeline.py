"""
Tests for the pipeline orchestrator and CLI.
All external steps are mocked so no real files or tools are needed.
"""

import yaml
import pytest
from pathlib import Path
from unittest.mock import patch, MagicMock
from click.testing import CliRunner

from grid.pipeline import run_wgs_pipeline
from grid.cli import cli

# ── helpers ────────────────────────────────────────────────────────────────


def write_config(tmp_path, extra=None):
    cfg = {
        "index": {"run": False},
        "count_reads": {"run": False},
        "mosdepth": {
            "run": False,
            "normalize": {"run": False},
            "neighbors": {"run": False},
            "gc_normalize": {"run": False},
        },
        "compute_diploid_genotypes": {"run": False},
        "compute_haploid_genotypes": {"run": False},
    }
    if extra:
        # Deep-merge mosdepth sub-keys so callers only need to override
        # what they care about without dropping gc_normalize.
        for key, val in extra.items():
            if key == "mosdepth" and isinstance(val, dict):
                cfg["mosdepth"].update(val)
            else:
                cfg[key] = val
    p = tmp_path / "config.yaml"
    p.write_text(yaml.dump(cfg))
    return str(p)


# ── pipeline tests ─────────────────────────────────────────────────────────


def test_pipeline_raises_without_config():
    with pytest.raises(Exception):
        run_wgs_pipeline(console=None, config=None)


def test_pipeline_all_steps_disabled(tmp_path):
    config = write_config(tmp_path)
    run_wgs_pipeline(console=None, config=config)


def test_pipeline_calls_check_index_when_index_disabled(tmp_path):
    config = write_config(tmp_path, {"index": {"run": False}})
    with patch("grid.utils.utils.check_index") as mock_ci:
        run_wgs_pipeline(console=None, config=config)
        mock_ci.assert_called_once()


def test_pipeline_calls_create_index_when_enabled(tmp_path):
    config = write_config(tmp_path, {"index": {"run": True}})
    with patch("grid.utils.utils.create_index") as mock_ci:
        run_wgs_pipeline(console=None, config=config)
        mock_ci.assert_called_once()


def test_pipeline_calls_count_reads_when_enabled(tmp_path):
    config = write_config(tmp_path, {"count_reads": {"run": True}})
    with patch("grid.utils.count_reads.count_reads") as mock_cr:
        run_wgs_pipeline(console=None, config=config)
        mock_cr.assert_called_once()


def test_pipeline_calls_mosdepth_when_enabled(tmp_path):
    config = write_config(tmp_path, {"mosdepth": {"run": True}})
    with patch("grid.utils.mosdepth.compute_mosdepth") as mock_md:
        run_wgs_pipeline(console=None, config=config)
        mock_md.assert_called_once()


def test_pipeline_calls_normalize_when_enabled(tmp_path):
    config = write_config(tmp_path, {"mosdepth": {"normalize": {"run": True}}})
    with patch("grid.utils.normalize_mosdepth.normalize_mosdepth") as mock_nm:
        run_wgs_pipeline(console=None, config=config)
        mock_nm.assert_called_once()


def test_pipeline_calls_find_neighbors_when_enabled(tmp_path):
    config = write_config(tmp_path, {"mosdepth": {"neighbors": {"run": True}}})
    with patch("grid.utils.find_neighbors.find_neighbors") as mock_fn:
        run_wgs_pipeline(console=None, config=config)
        mock_fn.assert_called_once()


def test_pipeline_calls_compute_diploid_when_enabled(tmp_path):
    config = write_config(tmp_path, {"compute_diploid_genotypes": {"run": True}})
    with patch("grid.utils.compute_dipcn.compute_diploid_genotypes") as mock_dip:
        run_wgs_pipeline(console=None, config=config)
        mock_dip.assert_called_once()


def test_pipeline_calls_hi_inference_when_enabled(tmp_path):
    config = write_config(tmp_path, {"compute_haploid_genotypes": {"run": True}})
    with patch("grid.utils.hi_inference.hi_inference") as mock_hi:
        run_wgs_pipeline(console=None, config=config)
        mock_hi.assert_called_once()


def test_pipeline_step_exception_does_not_crash_pipeline(tmp_path):
    config = write_config(tmp_path, {"count_reads": {"run": True}})
    with patch("grid.utils.count_reads.count_reads", side_effect=RuntimeError("boom")):
        run_wgs_pipeline(console=None, config=config)


# ── CLI tests ──────────────────────────────────────────────────────────────


def test_cli_wgs_command_runs(tmp_path):
    config = write_config(tmp_path)
    runner = CliRunner()
    with patch("grid.pipeline.run_wgs_pipeline") as mock_pipeline:
        result = runner.invoke(cli, ["wgs", config])
    assert result.exit_code == 0, result.output
    mock_pipeline.assert_called_once()


def test_cli_wgs_missing_config():
    runner = CliRunner()
    result = runner.invoke(cli, ["wgs", "/nonexistent/config.yaml"])
    assert result.exit_code != 0


def test_cli_help():
    runner = CliRunner()
    result = runner.invoke(cli, ["--help"])
    assert result.exit_code == 0
    assert "GRiD" in result.output


def test_cli_wgs_help():
    runner = CliRunner()
    result = runner.invoke(cli, ["wgs", "--help"])
    assert result.exit_code == 0


# ── compute_diploid_genotypes integration ──────────────────────────────────


def test_compute_diploid_genotypes_end_to_end(tmp_path):
    """Integration test: build minimal read-counts + neighbors files and
    run compute_diploid_genotypes.

    Output format (from source):
      - file: <output_file_prefix>.txt  (hardcoded .txt regardless of output_file_type)
      - content: space-separated, no header — one line per sample: "ID value"
      - neighbors file: <prefix>.zMax<zmax:g>.txt.gz  (capital M, :g format)
    """
    import gzip
    from contextlib import contextmanager
    from unittest.mock import patch, MagicMock
    from grid.utils.compute_dipcn import compute_diploid_genotypes

    # Write read counts file — no header, space-separated ID count
    counts_file = tmp_path / "counts.tsv"
    counts_file.write_text("S1 100\nS2 90\nS3 110\n")

    # Neighbors file: <prefix>.zMax2.txt.gz  (zmax=2, :g → "2")
    nbrs_file = tmp_path / "nbrs.zMax2.txt.gz"
    with gzip.open(nbrs_file, "wt") as f:
        f.write("S1\t1.0\tS2\t1.0\t0.1\tS3\t1.0\t0.2\n")
        f.write("S2\t1.0\tS1\t1.0\t0.1\tS3\t1.0\t0.15\n")
        f.write("S3\t1.0\tS1\t1.0\t0.2\tS2\t1.0\t0.15\n")

    config = {
        "output_dir": str(tmp_path),
        "output_file_type": "tsv",
        "compute_diploid_genotypes": {
            "output_file_prefix": "dipcn",
            "n_nbr": 10,
        },
        "count_reads": {"output_file_prefix": "counts"},
        "mosdepth": {
            "neighbors": {
                "zmax": 2,
                "output_file_prefix": "nbrs",
            },
            "gc_normalize": {"run": False},
        },
    }

    @contextmanager
    def _noop_progress(*args, **kwargs):
        mock_progress = MagicMock()
        mock_task = MagicMock()
        yield mock_progress, mock_task

    with patch("grid.utils.compute_dipcn.progress_bar", _noop_progress):
        compute_diploid_genotypes(config, console=None)

    # Output is always .txt, space-separated, no header
    out = tmp_path / "dipcn.txt"
    assert out.exists(), f"Expected {out} but not found. Files: {list(tmp_path.iterdir())}"
    lines = [l for l in out.read_text().strip().split("\n") if l.strip()]
    assert len(lines) == 3
    for line in lines:
        parts = line.split()
        assert len(parts) == 2
        assert parts[0] in ("S1", "S2", "S3")
        # value should be numeric or NA
        if parts[1] != "NA":
            float(parts[1])  # raises if not numeric
