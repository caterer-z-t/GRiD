import gzip
import os
import pytest
from pathlib import Path

from grid.utils.utils import (
    has_index,
    get_samples,
    find_file,
    setup_output_file,
    create_region_string,
    open_maybe_gz,
    get_flags,
)


# --- has_index ---

def test_has_index_cram_dot_crai(tmp_path):
    cram = tmp_path / "sample.cram"
    cram.touch()
    (tmp_path / "sample.cram.crai").touch()
    assert has_index(str(cram), "CRAM") is True

def test_has_index_cram_replaced_suffix(tmp_path):
    cram = tmp_path / "sample.cram"
    cram.touch()
    (tmp_path / "sample.crai").touch()
    assert has_index(str(cram), "CRAM") is True

def test_has_index_bam(tmp_path):
    bam = tmp_path / "sample.bam"
    bam.touch()
    (tmp_path / "sample.bam.bai").touch()
    assert has_index(str(bam), "BAM") is True

def test_has_index_missing(tmp_path):
    cram = tmp_path / "sample.cram"
    cram.touch()
    assert has_index(str(cram), "CRAM") is False

def test_has_index_invalid_type(tmp_path):
    f = tmp_path / "sample.xyz"
    f.touch()
    assert has_index(str(f), "XYZ") is False

def test_has_index_case_insensitive(tmp_path):
    cram = tmp_path / "sample.cram"
    cram.touch()
    (tmp_path / "sample.cram.crai").touch()
    assert has_index(str(cram), "cram") is True


# --- get_samples ---

def test_get_samples(tmp_path):
    samples_file = tmp_path / "samples.txt"
    samples_file.write_text("S1\nS2\nS3\n")
    assert get_samples(str(samples_file)) == ["S1", "S2", "S3"]

def test_get_samples_blank_lines(tmp_path):
    samples_file = tmp_path / "samples.txt"
    samples_file.write_text("S1\n\nS2\n\n")
    assert get_samples(str(samples_file)) == ["S1", "S2"]

def test_get_samples_empty(tmp_path):
    samples_file = tmp_path / "samples.txt"
    samples_file.write_text("")
    assert get_samples(str(samples_file)) == []


# --- find_file ---

def test_find_file_found(tmp_path):
    f = tmp_path / "mysample.cram"
    f.touch()
    result = find_file(str(tmp_path), "mysample", "cram")
    assert result is not None
    assert "mysample" in result

def test_find_file_not_found(tmp_path):
    result = find_file(str(tmp_path), "nosuchsample", "cram")
    assert result is None

def test_find_file_no_type(tmp_path):
    result = find_file(str(tmp_path), "sample", None)
    assert result is None


# --- setup_output_file ---

def test_setup_output_file(tmp_path):
    out = tmp_path / "out" / "result.tsv"
    path = setup_output_file(str(out), "chr6", 1000, 2000)
    assert path.exists()
    content = path.read_text()
    assert content == "Sample\tchr6:1000-2000\n"

def test_setup_output_file_creates_dirs(tmp_path):
    out = tmp_path / "deep" / "nested" / "result.tsv"
    setup_output_file(str(out), "chr1", 0, 500)
    assert out.exists()


# --- create_region_string ---

def test_create_region_string_valid():
    assert create_region_string("chr6", 100, 200) == "chr6:100-200"

def test_create_region_string_missing_chrom():
    with pytest.raises(ValueError):
        create_region_string(None, 100, 200)

def test_create_region_string_missing_start():
    with pytest.raises(ValueError):
        create_region_string("chr6", None, 200)


# --- open_maybe_gz ---

def test_open_maybe_gz_plain(tmp_path):
    f = tmp_path / "plain.txt"
    f.write_text("hello")
    with open_maybe_gz(str(f)) as fh:
        assert fh.read() == "hello"

def test_open_maybe_gz_gzipped(tmp_path):
    f = tmp_path / "data.txt.gz"
    with gzip.open(f, "wt") as fh:
        fh.write("hello gz")
    with open_maybe_gz(str(f)) as fh:
        assert fh.read() == "hello gz"


# --- get_flags ---

def test_get_flags():
    config = {"count_reads": {"flags": [83, 147, None]}}
    flags = get_flags(config, "count_reads")
    assert flags == [83, 147]

def test_get_flags_missing_key():
    assert get_flags({}, "count_reads") == []
