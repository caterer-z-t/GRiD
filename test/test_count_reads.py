# In[0]: Imports
from pathlib import Path
from GRiD.count_reads import count_reads

# In[1]: File paths
TEST_DIR = Path(__file__).parent / "files"
CRAM_FILE = TEST_DIR / "HG00152.hg38.cram"
REFERENCE_FASTA = TEST_DIR / "hg38.fa"

# In[2]: Tests
def test_count_reads_cram_region():
    """Test count_reads on a small region in the CRAM file."""
    region = "chr6:160605062-160605662"  # small region
    count = count_reads(str(CRAM_FILE), region, reference=str(REFERENCE_FASTA))
    
    # Expect some reads
    assert count > 0, f"Expected reads in region {region}, got {count}"