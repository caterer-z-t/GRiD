# LPA Realignment Module

Realign reads from the LPA KIV-2 VNTR region to a high-quality reference sequence for improved copy number estimation.

## Purpose

Standard genome alignment may be suboptimal in repeat regions due to:
- **Repetitive sequences** confusing aligners
- **Copy number variation** causing ambiguous mappings
- **Sequence homology** between repeat units

This module realigns reads specifically to the LPA KIV-2 reference, enabling:
1. **Precise read assignment** to specific exons (1A, 1B, KIV-3, etc.)
2. **Improved copy number accuracy** through targeted alignment
3. **Variant detection** within repeat units

## Main Module

### `align_lpa.py`
Located in parent directory: `grid/utils/align_lpa.py`

**Key Function:**
```python
def run_lpa_realignments(
    cram_dir: str,
    reference_fa: str,
    lpa_ref_fasta: str,
    positions_file: str,
    genome_build: str,
    chrom: str,
    start: int,
    end: int,
    output_file: str,
    threads: int = 1
) -> None
```

**Parameters:**
- `cram_dir`: Directory containing CRAM files
- `reference_fa`: Genome reference FASTA (for CRAM decoding)
- `lpa_ref_fasta`: LPA KIV-2 specific reference FASTA
- `positions_file`: TSV with exon positions in LPA reference
- `genome_build`: Genome build (hg19/hg37/hg38)
- `chrom`: Chromosome (e.g., "chr6")
- `start`: VNTR region start
- `end`: VNTR region end
- `output_file`: Output TSV with realignment counts
- `threads`: Parallel processing threads

**Output Format:**
```
sample_id       exon_1A exon_1B exon_1B_KIV3    exon_1B_notKIV3
sample001       234     1456    789             667
sample002       198     1523    823             700
sample003       267     1389    721             668
```

## Usage

### Via CLI
```bash
grid lpa-realign \
    --cram-dir data/CRAMs \
    --output-file output/lpa_realignment_counts.tsv \
    --ref-fasta refs/hs37d5.fa \
    --lpa-ref-fasta refs/lpa_kiv2_reference.fa \
    --positions-file refs/lpa_exon_positions.tsv \
    --genome-build hg38 \
    --chrom chr6 \
    --start 160500000 \
    --end 160600000 \
    --threads 4
```

### As Python Module
```python
from grid.utils.align_lpa import run_lpa_realignments

run_lpa_realignments(
    cram_dir="data/CRAMs",
    reference_fa="refs/hs37d5.fa",
    lpa_ref_fasta="refs/lpa_kiv2_reference.fa",
    positions_file="refs/lpa_exon_positions.tsv",
    genome_build="hg38",
    chrom="chr6",
    start=160500000,
    end=160600000,
    output_file="output/realignment_counts.tsv",
    threads=4
)
```

## LPA Reference Structure

The LPA gene contains multiple exons in the KIV-2 VNTR:

```
LPA Gene Structure:
[KIV-1][KIV-2 (variable copies)][KIV-3][KIV-4]...[Protease]

KIV-2 Repeats:
- Exon 1A: Diagnostic marker for KIV-2 copy number
- Exon 1B: Present in all KIV-2 copies
- Exon 1B_KIV3: Specific to KIV-3 (single copy)
- Exon 1B_notKIV3: Non-KIV3 copies
```

## Realignment Algorithm

### Step 1: Extract Reads
For each CRAM:
1. Fetch reads overlapping LPA VNTR region
2. Extract read sequences and qualities
3. Convert to FASTQ format

### Step 2: Realignment
Using BWA-MEM or similar aligner:
```bash
bwa mem lpa_reference.fa reads.fastq > aligned.sam
```

**Parameters:**
- Optimized for short, repetitive sequences
- Allowing multiple alignments per read
- Sensitive to detect all exon matches

### Step 3: Count Exon Hits
For each aligned read:
1. Check alignment position against exon coordinates
2. Assign read to exon type (1A, 1B, 1B_KIV3, etc.)
3. Count reads per exon type
4. Handle multi-mapping reads appropriately

### Step 4: Aggregate Results
- Combine counts across all samples
- Write to TSV output

## Positions File Format

The positions file defines exon boundaries in the LPA reference:

```
exon_type       start   end     description
exon_1A         1000    1150    KIV-2 specific marker
exon_1B         2000    2200    All KIV-2 copies
exon_1B_KIV3    3000    3200    KIV-3 specific
exon_1B_notKIV3 4000    4200    Non-KIV3 copies
```

**Columns:**
- `exon_type`: Exon identifier
- `start`: Start position in LPA reference (0-based)
- `end`: End position in LPA reference (0-based)
- `description`: Human-readable description

## Exon Types Explained

### Exon 1A
- **Uniqueness:** Found only in KIV-2 repeats
- **Copy number:** ~1 per KIV-2 repeat unit
- **Use:** Primary marker for KIV-2 CN estimation
- **Formula component:** Strong predictor (weight ~35 in final formula)

### Exon 1B
- **Uniqueness:** Present in all KIV copies (KIV-2, KIV-3, etc.)
- **Copy number:** Sum of all KIV types
- **Use:** Secondary marker
- **Formula component:** Moderate predictor (weight ~5)

### Exon 1B_KIV3
- **Uniqueness:** Specific to KIV-3 (single copy per haplotype)
- **Copy number:** 1-2 (diploid)
- **Use:** Normalization anchor point
- **Interpretation:** Should be ~2 in diploid individuals

### Exon 1B_notKIV3
- **Uniqueness:** KIV copies excluding KIV-3
- **Copy number:** Variable
- **Use:** Additional normalization
- **Calculation:** Exon_1B - Exon_1B_KIV3

## Copy Number Formula

The final KIV-2 copy number is estimated using:

```
diploid_CN = 34.9 × exon_1A + 5.2 × exon_1B - 1
haploid_CN = diploid_CN / 2
```

**Coefficients derived from:**
- Regression against orthogonal validation (PacBio long-reads, qPCR)
- Empirical calibration on reference samples
- Published in Mukamel et al. (2021)

## Technical Details

### Aligner Choice
**Recommended:** BWA-MEM
- Fast for short reads
- Handles repetitive sequences well
- Standard in genomics pipelines

**Alternative:** Bowtie2, minimap2
- May offer speed/sensitivity trade-offs

### Performance
- **Parallelization:** Process CRAMs in parallel
- **Speed:** ~2-5 minutes per CRAM
- **Memory:** <4GB per process
- **Disk:** Minimal (streaming)

### Dependencies
- **External:** `bwa`, `samtools`
- **Python:** `pysam`, `pandas`

## Quality Control

### Expected Counts (30X WGS)
For a sample with 20 KIV-2 copies (diploid):
```
exon_1A: ~600 reads (20 copies × ~30X)
exon_1B: ~750 reads (includes KIV-3)
exon_1B_KIV3: ~60 reads (2 copies)
```

### Validation Checks
1. **Exon 1B > Exon 1A:** Always true (includes KIV-3)
2. **Exon 1B_KIV3 ≈ 60 reads:** Should be ~2 diploid copies
3. **Read ratios consistent:** Within expected range

### Outlier Detection
**Warning signs:**
- Exon_1A = 0: Complete mapping failure
- Exon_1B_KIV3 > 200: Possible KIV-3 duplication
- Very low counts: Low coverage or alignment issues

## Troubleshooting

**No reads aligned:**
- Check LPA reference matches genome build
- Verify positions file coordinates
- Ensure BWA index exists for LPA reference

**Unexpected ratios:**
```python
# Check exon relationships
ratio = exon_1B / exon_1A
# Should be ~1.2-1.5 for most samples
```

**Index missing:**
```bash
# Build BWA index
bwa index lpa_reference.fa
```

**High KIV-3 counts:**
- May indicate KIV-3 duplication (rare variant)
- Verify with orthogonal methods
- Check family members if available

## References

- **Mukamel et al. (2021)** - Original LPA CN methodology
- **Hao et al. (2024)** - vntrwrap implementation
- **Li & Durbin (2009)** - BWA aligner algorithm

## Future Enhancements

- [ ] Support for long-read data (PacBio, ONT)
- [ ] Phased copy number estimation
- [ ] Variant calling within KIV-2 repeats
- [ ] Integration with graph-based alignment