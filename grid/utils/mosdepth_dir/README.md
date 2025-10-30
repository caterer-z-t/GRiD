# Mosdepth Coverage Analysis Module

Compute per-base or binned read depth coverage across genomic regions using the `mosdepth` tool.

## Purpose

This module generates coverage profiles for all samples, which are essential for:
1. Identifying high-coverage VNTR regions
2. Normalizing copy number estimates across samples
3. Quality control and outlier detection

## Main Module

### `run_mosdepth.py`
Located in parent directory: `grid/utils/run_mosdepth.py`

**Key Function:**
```python
def run_mosdepth(
    cram_dir: str,
    output_file: str,
    ref_fasta: str,
    chrom: str,
    start: int,
    end: int,
    region_name: str,
    work_dir: str = "~/mosdepth_work",
    by: int = 1000,
    fast_mode: bool = True,
    threads: int = 1
) -> None
```

**Parameters:**
- `cram_dir`: Directory with CRAM files
- `output_file`: Final aggregated TSV output
- `ref_fasta`: Reference genome FASTA
- `chrom`: Chromosome (e.g., "chr6")
- `start`: Region start position
- `end`: Region end position
- `region_name`: Name for mosdepth output (e.g., "LPA_VNTR")
- `work_dir`: Temporary directory for mosdepth outputs
- `by`: Bin size in base pairs (default: 1000)
- `fast_mode`: Use mosdepth --fast-mode (recommended: True)
- `threads`: Parallel processing threads

**Output Format:**
```
sample_id       position        depth
sample001       160000000       45.2
sample001       160001000       48.7
sample002       160000000       52.3
```

## Usage

### Via CLI
```bash
grid mosdepth \
    --cram-dir data/CRAMs \
    --output-file output/coverage.tsv \
    --ref-fasta refs/hs37d5.fa \
    --chrom chr6 \
    --start 160000000 \
    --end 160100000 \
    --region-name LPA_VNTR \
    --work-dir /tmp/mosdepth \
    --by 1000 \
    --fast \
    --threads 4
```

### As Python Module
```python
from grid.utils.run_mosdepth import run_mosdepth

run_mosdepth(
    cram_dir="data/CRAMs",
    output_file="output/coverage.tsv",
    ref_fasta="refs/hs37d5.fa",
    chrom="chr6",
    start=160000000,
    end=160100000,
    region_name="LPA_VNTR",
    by=1000,
    fast_mode=True,
    threads=4
)
```

## Mosdepth Integration

This module wraps the `mosdepth` tool with GRiD-specific configurations.

### Mosdepth Command Structure
```bash
mosdepth \
    --by <bin_size> \
    --chrom <chromosome> \
    --fast-mode \
    --threads <n> \
    <output_prefix> \
    <input.cram>
```

### Fast Mode
The `--fast-mode` flag is enabled by default:
- **Faster:** ~3-5x speed improvement
- **Trade-off:** May miss some edge cases in coverage calculation
- **Recommended:** Use for initial analysis and large datasets

## Output Files

### Intermediate Files (in work_dir)
For each sample, mosdepth generates:
```
<sample>.regions.bed.gz         # Per-region coverage summary
<sample>.mosdepth.global.dist.txt  # Global coverage distribution
<sample>.mosdepth.region.dist.txt  # Region coverage distribution
```

### Final Output
Aggregated TSV with columns:
- `sample_id`: Sample identifier (from CRAM filename)
- `position`: Genomic position (bin start)
- `depth`: Mean coverage depth in that bin

## Technical Details

### Algorithm
1. Create BED file for target region
2. For each CRAM in parallel:
   - Run mosdepth with specified parameters
   - Parse `.regions.bed.gz` output
   - Extract coverage values for target region
3. Aggregate all results into single TSV
4. Clean up intermediate files (optional)

### Performance
- **Parallelization:** Process multiple CRAMs simultaneously
- **Memory:** Low memory footprint with streaming
- **Speed:** ~1-2 minutes per CRAM for 100kb region
- **Disk:** Intermediate files are small (<1MB per sample)

### Dependencies
- **External:** `mosdepth` (must be in PATH)
- **Python:** `pandas`, `subprocess`, `pathlib`

## Binning Strategies

### Per-base (--by 0)
```bash
--by 0
```
- Most detailed coverage profile
- Large output files
- Slower processing
- Use for small regions (<10kb)

### Fixed bins (default)
```bash
--by 1000  # 1kb bins
```
- Balanced detail vs. speed
- Recommended for VNTR analysis
- Typical bin sizes: 500-2000bp

### Larger bins
```bash
--by 5000  # 5kb bins
```
- Faster processing
- Less detail
- Use for genome-wide analysis

## Quality Control

### Expected Coverage
For 30X WGS data in a VNTR region:
- Single-copy: ~30X coverage
- Duplicated: ~60X coverage
- High-copy (10+): 300X+ coverage

### Outlier Detection
Samples with unusual coverage patterns may indicate:
- Library preparation issues
- Alignment problems
- True structural variants

## Troubleshooting

**Mosdepth not found:**
```
Error: mosdepth: command not found
```
→ Install mosdepth: `conda install -c bioconda mosdepth`

**Empty output:**
- Check BED coordinates match reference
- Verify CRAM has coverage in region
- Ensure chromosome naming matches ("chr6" vs "6")

**Slow performance:**
- Use `--fast` mode
- Increase bin size (`--by`)
- Reduce `--threads` if I/O limited