# GRiD Utilities Module

Core computational utilities for the GRiD pipeline. Each module corresponds to a specific pipeline step and can be used independently or through the CLI.

## Directory Structure

```
utils/
├── align_lpa_dir/              # LPA realignment utilities
├── compute_dipcn_dir/          # Diploid copy number computation
├── count_reads_dir/            # Read counting utilities
├── estimate_kiv_dir/           # KIV-2 copy number estimation
├── find_neighbors_dir/         # Neighbor-finding algorithms
├── helper_dir/                 # Shared helper functions
├── mosdepth_dir/               # Coverage analysis utilities
├── normalize_mosdepth_dir/     # Coverage normalization
└── *.py                        # Main module scripts
```

## Core Modules

### File Management
- **`ensure_crai.py`** - Ensures CRAM index (.crai) files exist
- **`subset_cram.py`** - Subset CRAMs to specific genomic regions
- **`google_cloud_copy.py`** - (Deprecated) Google Cloud Storage utilities

### Read Counting & Coverage
- **`count_reads.py`** - Count properly paired reads in VNTR regions
- **`run_mosdepth.py`** - Compute per-base coverage using mosdepth
- **`run_normalized_mosdepth.py`** - Normalize coverage across samples

### Analysis Pipeline
- **`find_neighbors.py`** - Identify similar samples for normalization
- **`extract_reference.py`** - Extract VNTR reference sequences from genome
- **`align_lpa.py`** - Realign reads to LPA KIV-2 reference
- **`compute_dipcn.py`** - Calculate diploid copy numbers
- **`estimate_kiv.py`** - Estimate final KIV-2 copy numbers

## Usage Patterns

### As CLI Commands
All utilities are exposed through the `grid` CLI:
```bash
grid count-reads [options]
grid mosdepth [options]
grid normalize-mosdepth [options]
grid find-neighbors [options]
grid lpa-realign [options]
grid compute-dipcn [options]
grid estimate-kiv [options]
```

### As Python Modules
Utilities can also be imported directly:
```python
from grid.utils.run_mosdepth import run_mosdepth
from grid.utils.count_reads import count_reads

# Use programmatically
run_mosdepth(
    cram_dir="data/crams",
    output_file="out/coverage.tsv",
    ref_fasta="refs/hs37d5.fa",
    chrom="chr6",
    start=160000000,
    end=160100000
)
```

## Module Dependencies

Most modules require:
- `pysam` for BAM/CRAM handling
- `pandas` for data manipulation
- `numpy` for numerical operations
- External tools: `samtools`, `mosdepth`, `bwa`

See individual module directories for specific dependencies and detailed documentation.

## Helper Functions

Common utilities shared across modules are in `helper_dir/`:
- Region string formatting
- File discovery
- Result display
- Configuration loading
- Output file setup

See `helper_dir/README.md` for details.