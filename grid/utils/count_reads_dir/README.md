# Read Counting Module

Count properly paired reads within specified genomic regions (e.g., LPA KIV-2 VNTR) across multiple CRAM files.

## Purpose

This module is the first analytical step in the GRiD pipeline. It quantifies how many properly paired, mapped reads overlap with the repeat region of interest. These counts are essential for downstream copy number estimation.

## Main Module

### `count_reads.py`
Located in parent directory: `grid/utils/count_reads.py`

**Key Function:**
```python
def count_reads(
    cram_dir: str,
    output_file: str,
    ref_fasta: str,
    chrom: str,
    start: int,
    end: int,
    config: str,
    threads: int = 1
) -> None
```

**Parameters:**
- `cram_dir`: Directory containing CRAM files
- `output_file`: Output TSV file path
- `ref_fasta`: Reference genome FASTA (for CRAM decoding)
- `chrom`: Chromosome name (e.g., "chr6")
- `start`: Region start coordinate (0-based)
- `end`: Region end coordinate (0-based)
- `config`: YAML config with SAM flag filters
- `threads`: Parallel processing threads

**Output Format:**
```
sample_id       read_count
sample001       1234
sample002       2345
sample003       3456
```

## Usage

### Via CLI
```bash
grid count-reads \
    --cram-dir data/CRAMs \
    --output-file output/read_counts.tsv \
    --ref-fasta refs/hs37d5.fa \
    --chrom chr6 \
    --start 160000000 \
    --end 160100000 \
    --config config.yaml \
    --threads 4
```

### As Python Module
```python
from grid.utils.count_reads import count_reads

count_reads(
    cram_dir="data/CRAMs",
    output_file="output/counts.tsv",
    ref_fasta="refs/hs37d5.fa",
    chrom="chr6",
    start=160000000,
    end=160100000,
    config="config.yaml",
    threads=4
)
```

## SAM Flag Filtering

Reads are filtered based on flags specified in the YAML config:

```yaml
sam_flags:
  required_flags:
    - PROPER_PAIR    # Read is properly paired
    - READ_MAPPED    # Read is mapped
  excluded_flags:
    - DUPLICATE      # PCR/optical duplicate
    - SECONDARY      # Secondary alignment
    - SUPPLEMENTARY  # Supplementary alignment
```

**Common Required Flags:**
- `PROPER_PAIR` (0x2): Both reads in pair properly aligned
- `READ_MAPPED` (0x1): Read is mapped to reference

**Common Excluded Flags:**
- `DUPLICATE` (0x400): PCR or optical duplicate
- `SECONDARY` (0x100): Not primary alignment
- `SUPPLEMENTARY` (0x800): Supplementary alignment
- `QCFAIL` (0x200): Fails quality checks
- `UNMAPPED` (0x4): Read unmapped

## Technical Details

### Algorithm
1. Discover all CRAM files in specified directory
2. For each CRAM:
   - Ensure .crai index exists
   - Open CRAM with reference genome
   - Fetch reads overlapping target region
   - Filter reads by SAM flags
   - Count properly paired reads
3. Write results to TSV

### Performance
- **Parallelization:** Uses multiprocessing to process CRAMs concurrently
- **Memory:** Streams reads, minimal memory footprint
- **Speed:** ~30-60 seconds per CRAM (varies by region size and coverage)

### Dependencies
- `pysam` - CRAM/BAM reading
- `pandas` - Data output
- `multiprocessing` - Parallel execution

## Troubleshooting

**CRAM index missing:**
```
Error: index file <file>.crai does not exist
```
→ Run `grid crai` or ensure index files exist

**Reference mismatch:**
```
Error: [E::hts_open_format] fail to open file
```
→ Verify reference FASTA matches CRAM build (hg19/hg38)

**No reads found:**
- Check chromosome naming ("chr6" vs "6")
- Verify coordinates are correct for genome build
- Ensure region overlaps with sequencing data