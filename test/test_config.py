import pytest
from pathlib import Path

from grid.config import (
    validate_top_level,
    validate_steps,
    error_check_config,
    _get_nested,
    _is_enabled,
)


def minimal_config(tmp_path):
    samples = tmp_path / "samples.txt"
    samples.write_text("S1\n")
    ref = tmp_path / "ref.fa"
    ref.touch()
    return {
        "samples_file": str(samples),
        "directory_loc": str(tmp_path),
        "reference_genome": str(ref),
        "output_dir": str(tmp_path),
        "threads": 4,
        "file_type": "CRAM",
        "chrom": "chr6",
        "start_bp": 1000,
        "end_bp": 2000,
        "output_file_type": "tsv",
    }


# --- _get_nested ---

def test_get_nested_simple():
    assert _get_nested({"a": {"b": 1}}, "a", "b") == 1

def test_get_nested_missing():
    assert _get_nested({"a": {}}, "a", "b") is None

def test_get_nested_top_missing():
    assert _get_nested({}, "x") is None


# --- _is_enabled ---

def test_is_enabled_true():
    cfg = {"mosdepth": {"run": True}}
    assert _is_enabled(cfg, ("mosdepth",)) is True

def test_is_enabled_false():
    cfg = {"mosdepth": {"run": False}}
    assert _is_enabled(cfg, ("mosdepth",)) is False

def test_is_enabled_missing():
    assert _is_enabled({}, ("mosdepth",)) is False


# --- validate_top_level ---

def test_validate_top_level_ok(tmp_path):
    cfg = minimal_config(tmp_path)
    errors, warnings = [], []
    validate_top_level(cfg, errors, warnings)
    assert errors == []

def test_validate_top_level_missing_field(tmp_path):
    cfg = minimal_config(tmp_path)
    del cfg["chrom"]
    errors, warnings = [], []
    validate_top_level(cfg, errors, warnings)
    assert any("chrom" in e for e in errors)

def test_validate_top_level_wrong_type(tmp_path):
    cfg = minimal_config(tmp_path)
    cfg["threads"] = "four"
    errors, warnings = [], []
    validate_top_level(cfg, errors, warnings)
    assert any("threads" in e for e in errors)

def test_validate_top_level_missing_file(tmp_path):
    cfg = minimal_config(tmp_path)
    cfg["reference_genome"] = "/nonexistent/ref.fa"
    errors, warnings = [], []
    validate_top_level(cfg, errors, warnings)
    assert any("reference_genome" in e for e in errors)


# --- validate_steps (only runs checks for enabled steps) ---

def test_validate_steps_skips_disabled(tmp_path):
    cfg = minimal_config(tmp_path)
    # mosdepth.normalize not enabled → no error even if repeat_mask_file missing
    errors, warnings = [], []
    validate_steps(cfg, errors, warnings)
    assert errors == []

def test_validate_steps_required_field_missing(tmp_path):
    cfg = minimal_config(tmp_path)
    cfg["count_reads"] = {"run": True}  # flags is required but absent
    errors, warnings = [], []
    validate_steps(cfg, errors, warnings)
    assert any("flags" in e for e in errors)

def test_validate_steps_warns_on_default(tmp_path):
    cfg = minimal_config(tmp_path)
    cfg["count_reads"] = {"run": True, "flags": [83, 147]}
    # min_mapq not set → should warn
    errors, warnings = [], []
    validate_steps(cfg, errors, warnings)
    assert any("min_mapq" in w for w in warnings)


# --- error_check_config ---

def test_error_check_config_valid(tmp_path):
    cfg = minimal_config(tmp_path)
    error_check_config(cfg, console=None)  # should not raise

def test_error_check_config_raises_on_error(tmp_path):
    cfg = minimal_config(tmp_path)
    del cfg["threads"]
    with pytest.raises(ValueError, match="config error"):
        error_check_config(cfg, console=None)
