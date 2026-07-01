import gzip
import pytest
from pathlib import Path

from grid.utils.mosdepth import (
    compute_region_coverage,
    build_mosdepth_command,
    check_mosdepth_available,
    write_coverage_result,
    remove_intermediate_files,
)


# --- compute_region_coverage ---

def make_bed_gz(tmp_path, lines):
    f = tmp_path / "regions.bed.gz"
    with gzip.open(f, "wt") as fh:
        for line in lines:
            fh.write(line + "\n")
    return f

def test_compute_region_coverage_full_overlap(tmp_path):
    f = make_bed_gz(tmp_path, ["chr6\t1000\t2000\t30.0"])
    cov = compute_region_coverage(f, "chr6", 1000, 2000)
    assert cov == int(round(100 * 30.0))

def test_compute_region_coverage_partial_overlap(tmp_path):
    # bin covers 0–2000, region is 1000–2000 → 1000 bp overlap, depth=40
    f = make_bed_gz(tmp_path, ["chr6\t0\t2000\t40.0"])
    cov = compute_region_coverage(f, "chr6", 1000, 2000)
    assert cov == int(round(100 * 40.0))

def test_compute_region_coverage_multiple_bins(tmp_path):
    f = make_bed_gz(tmp_path, [
        "chr6\t1000\t2000\t20.0",
        "chr6\t2000\t3000\t40.0",
    ])
    cov = compute_region_coverage(f, "chr6", 1000, 3000)
    # weighted mean: (20*1000 + 40*1000) / 2000 = 30
    assert cov == int(round(100 * 30.0))

def test_compute_region_coverage_wrong_chrom(tmp_path):
    f = make_bed_gz(tmp_path, ["chr1\t1000\t2000\t30.0"])
    cov = compute_region_coverage(f, "chr6", 1000, 2000)
    assert cov == 0

def test_compute_region_coverage_no_overlap(tmp_path):
    f = make_bed_gz(tmp_path, ["chr6\t5000\t6000\t30.0"])
    cov = compute_region_coverage(f, "chr6", 1000, 2000)
    assert cov == 0

def test_compute_region_coverage_empty_file(tmp_path):
    f = make_bed_gz(tmp_path, [])
    cov = compute_region_coverage(f, "chr6", 1000, 2000)
    assert cov == 0


# --- build_mosdepth_command ---

def test_build_mosdepth_command_basic(tmp_path):
    cmd = build_mosdepth_command("sample.cram", "ref.fa", tmp_path / "out", 1000, False, 1)
    assert cmd[0] == "mosdepth"
    assert "--by" in cmd
    assert "1000" in cmd
    assert "-f" in cmd
    assert "--fast-mode" not in cmd

def test_build_mosdepth_command_fast_mode(tmp_path):
    cmd = build_mosdepth_command("sample.cram", "ref.fa", tmp_path / "out", 1000, True, 4)
    assert "--fast-mode" in cmd
    assert cmd[0] == "mosdepth"

def test_build_mosdepth_command_threads(tmp_path):
    cmd = build_mosdepth_command("s.cram", "r.fa", tmp_path / "o", 500, False, 8)
    idx = cmd.index("-t")
    assert cmd[idx + 1] == "8"


# --- write_coverage_result ---

def test_write_coverage_result(tmp_path):
    import threading
    f = tmp_path / "out.tsv"
    f.write_text("Sample\tchr6:0-1000\n")
    lock = threading.Lock()
    write_coverage_result(f, "S1", 3000, lock)
    content = f.read_text()
    assert "S1\t3000" in content


# --- remove_intermediate_files ---

def test_remove_intermediate_files(tmp_path):
    (tmp_path / "sample.mosdepth.global.dist.txt").touch()
    (tmp_path / "sample.regions.bed.gz.csi").touch()
    (tmp_path / "keep_me.tsv").touch()
    remove_intermediate_files(tmp_path)
    assert not (tmp_path / "sample.mosdepth.global.dist.txt").exists()
    assert not (tmp_path / "sample.regions.bed.gz.csi").exists()
    assert (tmp_path / "keep_me.tsv").exists()

def test_remove_intermediate_files_include_bed_gz(tmp_path):
    (tmp_path / "sample.regions.bed.gz").touch()
    remove_intermediate_files(tmp_path, include_region_bed_gz=True)
    assert not (tmp_path / "sample.regions.bed.gz").exists()

def test_remove_intermediate_files_empty_dir(tmp_path):
    remove_intermediate_files(tmp_path)  # should not raise


# --- check_mosdepth_available ---

def test_check_mosdepth_available_missing(monkeypatch):
    import shutil
    monkeypatch.setattr(shutil, "which", lambda _: None)
    with pytest.raises(RuntimeError, match="mosdepth"):
        check_mosdepth_available()

def test_check_mosdepth_available_present(monkeypatch):
    import shutil
    monkeypatch.setattr(shutil, "which", lambda _: "/usr/bin/mosdepth")
    check_mosdepth_available()  # should not raise
