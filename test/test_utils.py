# In[0]: Imports
import pytest
import shutil
from pathlib import Path
from GRiD.utils.utils import format_region, arg_check, ensure_index
import argparse

# In[1]: File paths
TEST_DIR = Path(__file__).parent / "files"
CRAM_FILE = TEST_DIR / "HG00152.hg38.cram"
CRAM_INDEX = TEST_DIR / "HG00152.hg38.cram.crai"
REFERENCE_FASTA = TEST_DIR / "hg38.fa"

# In[2]: Tests
def test_format_region_adds_chr():
    # Test chromosome formatting
    args = argparse.Namespace(chromosome="6", start=100, end=200)
    region = format_region(args)
    assert region == "chr6:100-200"

def test_format_region_keeps_chr():
    # Test chromosome already starts with chr
    args = argparse.Namespace(chromosome="chrX", start=500, end=600)
    region = format_region(args)
    assert region == "chrX:500-600"

def test_arg_check_missing_input(tmp_path):
    # Input file does not exist
    args = argparse.Namespace(input=str(tmp_path / "fake.bam"), region=None,
                              start=100, end=200, chromosome="1", reference=None)
    with pytest.raises(FileNotFoundError):
        arg_check(args)

def test_arg_check_cram_without_reference(tmp_path):
    # CRAM file without reference should raise ValueError
    fake_cram = tmp_path / "sample.cram"
    fake_cram.touch()
    args = argparse.Namespace(input=str(fake_cram), region=None,
                              start=100, end=200, chromosome="1", reference=None)
    with pytest.raises(ValueError):
        arg_check(args)

def test_arg_check_bam_index_missing(tmp_path):
    # BAM file without .bai should raise FileNotFoundError
    fake_bam = tmp_path / "sample.bam"
    fake_bam.touch()
    args = argparse.Namespace(input=str(fake_bam), region=None,
                              start=100, end=200, chromosome="1", reference=None)
    with pytest.raises(FileNotFoundError):
        arg_check(args)

def test_arg_check_region_missing(tmp_path):
    # Create a fake BAM file so the file check passes
    fake_bam = tmp_path / "sample.bam"
    fake_bam_index = tmp_path / "sample.bam.bai"
    fake_bam.touch()
    fake_bam_index.touch()

    # Missing region and missing chromosome/start/end
    args = argparse.Namespace(
        input=str(fake_bam),
        region=None,
        start=None,
        end=None,
        chromosome=None,
        reference=None
    )

    with pytest.raises(ValueError):
        arg_check(args)

def test_format_region_real():
    """Test region formatting from chromosome, start, end."""
    args = argparse.Namespace(chromosome="6", start=160605062, end=160605662)
    region = format_region(args)
    assert region == "chr6:160605062-160605662"

def test_arg_check_real_files():
    """Test arg_check does not raise errors for real files."""
    args = argparse.Namespace(
        input=str(CRAM_FILE),
        region="chr6:160605062-160605662",
        start=None,
        end=None,
        chromosome=None,
        reference=str(REFERENCE_FASTA),
    )
    # Should not raise any exception
    arg_check(args)

# In[3]: Tests for ensure_index
def test_ensure_index_file_not_found():
    """Raises FileNotFoundError if input file does not exist."""
    with pytest.raises(FileNotFoundError):
        ensure_index("nonexistent.cram")


def test_ensure_index_invalid_extension(tmp_path):
    """Raises ValueError if input file is not .bam or .cram."""
    txt_file = tmp_path / "fake.txt"
    txt_file.write_text("not a cram")
    with pytest.raises(ValueError):
        ensure_index(str(txt_file))