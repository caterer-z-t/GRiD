# Mosdepth Coverage Normalization Module

Normalize raw coverage data across samples to account for sequencing depth variation and technical artifacts.

## Purpose

Raw coverage values vary widely between samples due to:
- Different sequencing depths (20X vs 60X)
- Library preparation batch effects
- GC content bias
- Mappability differences

This module normalizes coverage using z-score transformation to enable fair comparison across samples.

## Main Module

### `run_normalized_mosdepth.py`
Located in parent directory: `grid/utils/run_normalized_mosdepth.py`

**Key Function:**
```python
def run_normalize_mosdepth(
    mosdepth_dir: str,
    output_file: str,
    repeat_mask: str,
    chrom: str,
    start: int,
    end: int,
    min_depth: int = 20,
    max_depth: int = 100,
    top_frac: float = 0.1,
    threads: int = 1
) -> None
```

**Parameters:**
- `mosdepth_dir`: Directory with mosdepth `.regions.bed.gz` files
- `output_file`: Output normalized z-depth file (`.gz`)
- `repeat_mask`: BED file of repeat-masked regions to exclude
- `chrom`: Chromosome
- `start`: Region start
- `end`: Region end
- `min_depth`: Minimum mean depth threshold (default: 20)
- `max_depth`: Maximum mean depth threshold (default: 100)
- `top_frac`: Top fraction of high-variance regions to keep (default: 0.1)
- `threads`: Parallel processing threads

**Output Format:**
Gzipped TSV with z-transformed depths:
```
sample_id       position        z_depth
sample001       160000000       1.23
sample001       160001000       -0.45
sample002       160000000       0.87
```

## Usage

### Via CLI
```bash
grid normalize-mosdepth \
    --mosdepth-dir output/mosdepth/ \
    --output-file output/normalized.zdepth.gz \
    --repeat-mask refs/repeat_mask.bed \
    --chrom chr6 \
    --start 160000000 \
    --end 160100000 \
    --min-depth 20 \
    --max-depth 100 \
    --top-frac 0.1 \
    --threads 4
```

### As Python Module
```python
from grid.utils.run_normalized_mosdepth import run_normalize_mosdepth

run_normalize_mosdepth(
    mosdepth_dir="output/mosdepth/",
    output_file="output/normalized.zdepth.gz",
    repeat_mask="refs/repeat_mask.bed",
    chrom="chr6",
    start=160000000,
    end=160100000,
    min_depth=20,
    max_depth=100,
    top_frac=0.1,
    threads=4
)
```

## Normalization Algorithm

### Step 1: Load Coverage Data
- Read all mosdepth `.regions.bed.gz` files
- Extract coverage for target region
- Create sample × position matrix

### Step 2: Filter Regions
**Repeat Masking:**
- Remove positions overlapping repeat-masked regions
- Prevents bias from low-complexity sequences

**Depth Filtering:**
- Remove positions with mean depth < `min_depth` (low coverage)
- Remove positions with mean depth > `max_depth` (extreme outliers)

**Variance Selection:**
- Calculate variance across samples for each position
- Keep top `top_frac` of positions by variance
- Focuses on informative, variable regions

### Step 3: Z-Score Transformation
For each sample independently:

```
z_depth = (depth - mean_depth) / std_depth
```

Where:
- `depth`: Raw coverage at position
- `mean_depth`: Sample's mean coverage across all positions
- `std_depth`: Sample's standard deviation

**Result:**
- Mean = 0, StdDev = 1 for each sample
- Comparable across samples with different depths

## Output Interpretation

### Z-Depth Values
- **z = 0**: Average coverage for this sample
- **z = +2**: 2 standard deviations above average (potential duplication)
- **z = -2**: 2 standard deviations below average (potential deletion)
- **z > +3**: Very high confidence variant

### Example Interpretation
```
Sample    Position      Z-Depth    Interpretation
S001      160050000     2.5        Likely duplication
S001      160051000     2.8        Strong duplication signal
S001      160052000     -0.3       Normal coverage
S002      160050000     0.1        Normal coverage
```

Sample S001 likely has a duplication in the 160050000-160051000 region.

## Parameter Tuning

### `min_depth` (default: 20)
- **Purpose:** Filter low-coverage positions
- **Too low:** Include unreliable positions
- **Too high:** Remove informative data
- **Recommendation:** 
  - 30X WGS: use 15-20
  - 60X WGS: use 30-40

### `max_depth` (default: 100)
- **Purpose:** Remove extreme outliers
- **Too low:** Filter true high-copy variants
- **Too high:** Include technical artifacts
- **Recommendation:**
  - Depends on expected max copy number
  - LPA KIV-2: 100-150 (max ~50 copies × 2 haploids)

### `top_frac` (default: 0.1)
- **Purpose:** Focus on most variable positions
- **Range:** 0.05 - 0.3
- **Lower values:** More stringent, fewer positions
- **Higher values:** More positions, may include noise
- **Recommendation:**
  - VNTR analysis: 0.1 (10% most variable)
  - CNV discovery: 0.2 (20% most variable)

## Repeat Mask BED Format

The repeat mask file should be a standard BED file:
```
chr6    160000000    160001000    SimpleRepeat
chr6    160005000    160005500    LINE
chr6    160010000    160011000    Satellite
```

**Purpose:**
- Exclude low-complexity sequences
- Remove problematic alignment regions
- Focus on informative VNTR sequences

**How to generate:**
1. Download RepeatMasker track from UCSC Genome Browser
2. Filter for target region
3. Optionally keep only certain repeat types

## Technical Details

### Performance
- **Parallelization:** Sample-level parallel processing
- **Memory:** Entire dataset loaded into memory (~GB for 1000 samples)
- **Speed:** ~5-10 minutes for 1000 samples, 100kb region

### Dependencies
- `pandas`, `numpy` - Data processing
- `scipy` - Statistical functions
- `multiprocessing` - Parallel execution

## Quality Control

### Expected Results
After normalization:
- All samples have mean ≈ 0, std ≈ 1
- High-copy regions show consistently high z-scores
- Low-noise regions show z-scores near 0

### Diagnostic Plots (manual)
```python
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("normalized.zdepth.gz", sep="\t")
sample_means = df.groupby("sample_id")["z_depth"].mean()

plt.hist(sample_means, bins=50)
plt.xlabel("Mean Z-Depth")
plt.ylabel("Samples")
plt.title("Should be centered at 0")
plt.show()
```

## Troubleshooting

**All z-scores near zero:**
- Region may not have sufficient variation
- Try reducing `top_frac` threshold
- Check that region contains actual VNTRs

**Extreme z-scores (|z| > 10):**
- May indicate technical artifacts
- Check for alignment issues
- Verify repeat mask is correct

**Memory errors:**
- Reduce number of samples processed at once
- Increase system RAM
- Use smaller genomic region