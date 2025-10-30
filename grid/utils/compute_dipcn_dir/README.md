# Diploid Copy Number Computation Module

Compute diploid copy numbers for LPA KIV-2 exons using neighbor-based normalization. This is the penultimate step before final KIV-2 CN estimation.

## Purpose

This module transforms raw realignment read counts into normalized diploid copy number estimates by:

1. **Neighbor-based normalization** - Compare each sample to similar individuals
2. **Batch effect correction** - Remove technical variation
3. **Diploid CN estimation** - Account for both haplotypes
4. **Exon-specific quantification** - Separate estimates for 1A, 1B, and subtypes

The output provides refined CN estimates for each exon type, which are combined in the final step to estimate total KIV-2 copy number.

## Main Module

### `compute_dipcn.py`
Located in parent directory: `grid/utils/compute_dipcn.py`

**Key Function:**
```python
def compute_dipcn_pipeline(
    count_file: Path,
    neighbor_file: Path,
    output_prefix: str,
    n_neighbors: int = 200,
    console: Optional[Console] = None
) -> None
```

**Parameters:**
- `count_file`: Realignment output TSV (from `lpa-realign`)
- `neighbor_file`: Neighbor file (from `find-neighbors`)
- `output_prefix`: Prefix for output files (e.g., "output/sample")
- `n_neighbors`: Number of top neighbors to use (default: 200)
- `console`: Rich console for pretty output (optional)

**Output Files:**
Four files are generated with diploid CN for each exon type:
```
<output_prefix>.exon1A.dipCN.txt
<output_prefix>.exon1B.dipCN.txt
<output_prefix>.exon1B_KIV3.dipCN.txt
<output_prefix>.exon1B_notKIV3.dipCN.txt
```

**Output Format (each file):**
```
sample_id       diploid_CN
sample001       18.45
sample002       22.13
sample003       15.89
```

## Usage

### Via CLI
```bash
grid compute-dipcn \
    --count-file output/lpa_realignment_counts.tsv \
    --neighbor-file output/neighbors.zMax2.0.txt.gz \
    --output-prefix output/sample_cohort \
    --n-neighbors 200
```

**Output:**
```
output/sample_cohort.exon1A.dipCN.txt
output/sample_cohort.exon1B.dipCN.txt
output/sample_cohort.exon1B_KIV3.dipCN.txt
output/sample_cohort.exon1B_notKIV3.dipCN.txt
```

### As Python Module
```python
from pathlib import Path
from grid.utils.compute_dipcn import compute_dipcn_pipeline

compute_dipcn_pipeline(
    count_file=Path("output/lpa_realignment_counts.tsv"),
    neighbor_file=Path("output/neighbors.zMax2.0.txt.gz"),
    output_prefix="output/sample_cohort",
    n_neighbors=200
)
```

## Algorithm Details

### Step 1: Load Data
1. **Read counts:** Load realignment counts for all samples
2. **Neighbors:** Load neighbor lists for each sample
3. **Validate:** Ensure all samples have neighbor information

### Step 2: Neighbor-Based Normalization

For each sample and each exon type:

```python
# Get top N neighbors
neighbors = neighbor_list[:n_neighbors]

# Calculate neighbor statistics
neighbor_counts = counts[neighbors]
neighbor_median = median(neighbor_counts)
neighbor_std = std(neighbor_counts)

# Normalize sample count
normalized_count = sample_count / neighbor_median
```

**Why neighbor-based?**
- **Reduces batch effects:** Similar samples likely from same batch
- **Population-aware:** Accounts for ancestry-related variation
- **Adaptive:** Each sample gets custom reference baseline

### Step 3: Copy Number Estimation

Convert normalized counts to diploid CN:

```python
# Assuming KIV-3 is diploid (2 copies) as anchor
reference_diploid = 2.0
diploid_CN = normalized_count r'\times{×}' reference_diploid
```

**For exon_1B_KIV3:**
- Expected to be 2 copies in diploid individuals (1 per haplotype)
- Used as normalization anchor
- Deviations may indicate rare KIV-3 variants

**For exon_1A and exon_1B:**
- Variable copy number
- Normalized against neighbors with similar genomic background
- Diploid estimate accounts for both haplotypes

### Step 4: Quality Control

Filter and flag potential issues:
- Extreme outliers (|z-score| > 4)
- Very low read counts (< 10 reads)
- Samples with insufficient neighbor data

## Copy Number Interpretation

### Exon 1A Diploid CN
**Range:** Typically 10-40 (diploid)
- **10-15:** Lower copy number genotype
- **20-30:** Common range
- **35-40:** Higher copy number genotype
- **>40:** Rare, verify with orthogonal methods

**Biological meaning:**
Each KIV-2 repeat unit has 1 exon 1A, so:
```
haploid KIV-2 copies ≈ diploid_CN / 2
```

### Exon 1B Diploid CN
**Range:** Typically 12-45 (diploid)
- Includes all KIV types (KIV-2, KIV-3, etc.)
- Should be slightly higher than exon 1A
- Difference reflects non-KIV-2 copies

**Relationship:**
```
exon_1B_CN ≈ exon_1A_CN + exon_1B_KIV3_CN
```

### Exon 1B_KIV3 Diploid CN
**Expected:** ~2.0 (diploid, single copy per haplotype)
- **1.5-2.5:** Normal range with measurement noise
- **<1.0:** Possible deletion (very rare)
- **>3.0:** Possible duplication (very rare)

**Use as QC:**
Samples with very aberrant KIV-3 CN may have:
- Alignment issues
- True rare variants
- Technical artifacts

### Exon 1B_notKIV3 Diploid CN
**Calculation:**
```
exon_1B_notKIV3 = exon_1B - exon_1B_KIV3
```

**Range:** Should closely match exon 1A
- Additional normalization check
- Validates exon relationships

## Parameter Tuning

### `n_neighbors` (default: 200)
**Purpose:** Number of similar samples for normalization baseline

- **Fewer (50-100):**
  - More stringent
  - Only closest matches
  - Better for homogeneous cohorts
  
- **More (300-500):**
  - More robust to outliers
  - Better for diverse cohorts
  - May include less similar samples

**Recommendation:**
- Small cohort (<500): 100-150
- Medium cohort (500-2000): 200-300
- Large cohort (>2000): 200-500

**Trade-off:**
- Too few: Sensitive to outlier neighbors
- Too many: Include dissimilar samples, reduce specificity

## Expected Results

### Typical Distribution
For 1000 samples from diverse populations:

```
Exon 1A Diploid CN:
  Mean: 25.3
  Median: 24.8
  Range: 12.5 - 42.7
  
Exon 1B_KIV3 Diploid CN:
  Mean: 2.04
  Median: 2.01
  Range: 1.7 - 2.3
```

### Quality Metrics

**Coefficient of Variation (CV):**
```python
cv = std(diploid_CN) / mean(diploid_CN)
```

- **Good:** CV < 0.3 for exon 1A
- **Acceptable:** CV < 0.5
- **Poor:** CV > 0.5 (check for batch effects)

**Exon Correlation:**
```python
correlation(exon_1A, exon_1B)
```

- **Good:** r > 0.85
- **Acceptable:** r > 0.75
- **Poor:** r < 0.75 (alignment or normalization issues)

## Validation

### Internal Consistency Checks

**Check 1: Exon relationships**
```python
# exon_1B should be >= exon_1A
assert (exon_1B >= exon_1A).all()

# exon_1B_KIV3 should be ~2
assert 1.5 < exon_1B_KIV3.mean() < 2.5
```

**Check 2: Correlation**
```python
import pandas as pd

df = pd.DataFrame({
    'exon_1A': exon_1A_CN,
    'exon_1B': exon_1B_CN
})

correlation = df.corr()
print(correlation)
# Should show strong positive correlation
```

### Comparison with Known Samples

If validation samples with known CN are available:
```python
known_samples = {
    'sample001': 20,  # Known KIV-2 CN
    'sample002': 25,
    'sample003': 15
}

for sample, true_CN in known_samples.items():
    estimated_CN = diploid_CN_df.loc[sample, 'diploid_CN'] / 2
    error = abs(estimated_CN - true_CN)
    print(f"{sample}: True={true_CN}, Est={estimated_CN:.1f}, Error={error:.1f}")
```

## Technical Details

### Performance
- **Memory:** <4GB for 1000 samples
- **Speed:** ~30-60 seconds for 1000 samples
- **Parallelization:** Not currently implemented (fast enough)

### Dependencies
- `pandas` - Data manipulation
- `numpy` - Numerical operations
- `scipy` - Statistical functions (optional)

## Troubleshooting

### All diploid CNs near 2
**Possible causes:**
- Normalization issue - check neighbor file
- All samples similar - low biological variation
- Incorrect reference scaling

**Solution:**
- Verify neighbor file is correct
- Check exon count distributions
- Ensure exon_1B_KIV3 used as anchor

### Extreme outliers (CN > 100)
**Possible causes:**
- Technical artifacts
- Very poor neighbors for that sample
- True extreme variant (rare)

**Solution:**
- Check read counts for that sample
- Inspect neighbor similarity
- Consider removing sample from analysis

### Negative diploid CNs
**Should never happen** - indicates serious error
- Check input count file format
- Verify neighbor file integrity
- Report as bug if reproducible

### Poor correlation between exons
**Possible causes:**
- Alignment issues in realignment step
- Batch effects not properly normalized
- Incorrect exon positions file

**Solution:**
- Rerun realignment with validated positions file
- Check for batch confounding in neighbor selection
- Increase `n_neighbors` for more robust normalization

## Next Step

After computing diploid CNs, proceed to final KIV-2 estimation:

```bash
grid estimate-kiv \
    --exon1a output/sample_cohort.exon1A.dipCN.txt \
    --exon1b output/sample_cohort.exon1B.dipCN.txt \
    --output output/final_kiv2_estimates.tsv
```

This combines the exon diploid CNs using the empirically-derived formula to produce final KIV-2 copy number estimates.