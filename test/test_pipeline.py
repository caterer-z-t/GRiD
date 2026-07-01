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
        "index":                    {"run": False},
        "count_reads":              {"run": False},
        "mosdepth":                 {"run": False, "normalize": {"run": False}, "neighbors": {"run": False}},
        "compute_diploid_genotypes": {"run": False},
        "compute_haploid_genotypes": {"run": False},
    }
    if extra:
        cfg.update(extra)
    p = tmp_path / "config.yaml"
    p.write_text(yaml.dump(cfg))
    return str(p)


# ── pipeline tests ─────────────────────────────────────────────────────────

def test_pipeline_raises_without_config():
    with pytest.raises(Exception):
        run_wgs_pipeline(console=None, config=None)

def test_pipeline_all_steps_disabled(tmp_path):
    config = write_config(tmp_path)
    # Should complete without error when every step is disabled
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
    config = write_config(tmp_path, {"mosdepth": {"run": True, "normalize": {"run": False}, "neighbors": {"run": False}}})
    with patch("grid.utils.mosdepth.compute_mosdepth") as mock_md:
        run_wgs_pipeline(console=None, config=config)
        mock_md.assert_called_once()

def test_pipeline_calls_normalize_when_enabled(tmp_path):
    config = write_config(tmp_path, {
        "mosdepth": {"run": False, "normalize": {"run": True}, "neighbors": {"run": False}}
    })
    with patch("grid.utils.normalize_mosdepth.normalize_mosdepth") as mock_nm:
        run_wgs_pipeline(console=None, config=config)
        mock_nm.assert_called_once()

def test_pipeline_calls_find_neighbors_when_enabled(tmp_path):
    config = write_config(tmp_path, {
        "mosdepth": {"run": False, "normalize": {"run": False}, "neighbors": {"run": True}}
    })
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
        # Pipeline catches exceptions per step — should not raise
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
