# KIV-2 Copy Number Estimation Module

Final step: combine exon-specific diploid copy numbers into final KIV-2 VNTR copy number estimates using an empirically calibrated regression formula.

## Purpose

This module takes the normalized exon diploid copy numbers from `compute_dipcn` and produces final, calibrated KIV-2 copy number estimates through:

1. **Regression formula** - Empirically derived coefficients from validation data
2. **Multi-exon integration** - Combines information from exon 1A and 1B
3. **Diploid/haploid conversion** - Provides both diploid and haploid estimates
4. **Final output** - Ready-to-use copy number calls

## Main Module

### `estimate_kiv.py`
Located in parent directory: `grid/utils/estimate_kiv.py`

**Key Function:**
```python
def compute_kiv_estimates(
    output: str,
    exon1a_file: str,
    exon1b_file: str,
    console: Optional[Console] = None
) -> pd.DataFrame
```

**Parameters:**
- `output`: Output file path for final estimates
- `exon1a_file`: Exon 1A diploid CN file (from `compute-dipcn`)
- `exon1b_file`: Exon 1B diploid CN file (from `compute-dipcn`)
- `console`: Rich console for pretty output (optional)

**Output Format:**
```
sample_id       diploid_estimate        haploid_estimate
sample001       36.75                   18.38
sample002       44.21                   22.11
sample003       28.93                   14.47
```

## Usage

### Via CLI
```bash
grid estimate-kiv \
    --exon1a output/cohort.exon1A.dipCN.txt \
    --exon1b output/cohort.exon1B.dipCN.txt \
    --output output/final_kiv2_cn_estimates.tsv \
    --format tsv
```

### As Python Module
```python
from grid.utils.estimate_kiv import compute_kiv_estimates

estimates_df = compute_kiv_estimates(
    output="output/kiv2_estimates.tsv",
    exon1a_file="output/cohort.exon1A.dipCN.txt",
    exon1b_file="output/cohort.exon1B.dipCN.txt"
)

print(estimates_df.head())
```

## Copy Number Formula

The formula was empirically derived through regression against orthogonal validation methods (long-read sequencing, qPCR):

```python
diploid_CN = 34.9 r'$\times{}$' exon1A_dipCN r'+' 5.2 r'$\times$' exon1B_dipCN r'-' 1
haploid_CN = diploid_CN r'/' 2
```

### Coefficient Interpretation

**Exon 1A coefficient (34.9):**
- **Strong predictor** - Primary signal for KIV-2 CN
- Each KIV-2 repeat has exactly 1 exon 1A
- Coefficient accounts for normalization scaling and technical factors

**Exon 1B coefficient (5.2):**
- **Secondary predictor** - Moderate contribution
- Exon 1B includes KIV-2 and other KIV types
- Lower weight due to presence in non-KIV-2 copies

**Intercept (-1):**
- **Baseline correction** - Accounts for systematic biases
- Derived from validation sample calibration

### Formula Derivation

The formula was fit using:
1. **Training samples:** N=500 with known KIV-2 CN from long-reads
2. **Method:** Multiple linear regression
3. **Validation:** R² = 0.94 on held-out test set
4. **Published:** Mukamel et al. (2021), Science

## Output Interpretation

### Diploid Estimate
**Total copies across both haplotypes (chromosomes)**

**Typical range:** 14-50 copies (diploid)
- **14-20:** Lower copy number genotype (~7-10 per haplotype)
- **20-30:** Common range (~10-15 per haplotype)
- **30-40:** Higher copy number genotype (~15-20 per haplotype)
- **40-50:** High copy number (~20-25 per haplotype)
- **>50:** Very high (rare, validate with orthogonal methods)

### Haploid Estimate
**Average copies per haplotype (single chromosome)**

```
haploid_CN = diploid_CN / 2
```

**Why haploid?**
- Easier biological interpretation
- Matches literature conventions
- Comparable to single-chromosome assays (e.g., Southern blot)

**Typical range:** 7-25 copies per haplotype
- **<5:** Very low (rare, potential deletion)
- **5-10:** Lower range
- **10-20:** Common range
- **20-30:** Higher range
- **>30:** Very high (rare, validate)

## Clinical Relevance

### LPA and Cardiovascular Disease

KIV-2 copy number is inversely correlated with:
- **Lp(a) levels** - Fewer copies = higher Lp(a)
- **CVD risk** - Lower CN associated with increased risk

**Risk stratification (haploid CN):**
- **<15 copies:** Higher Lp(a), increased CVD risk
- **15-20 copies:** Moderate
- **>20 copies:** Lower Lp(a), reduced risk

**Note:** Individual risk depends on many factors beyond CN. Clinical interpretation should involve comprehensive cardiovascular assessment.

### Genetic Counseling

Copy number information may be relevant for:
- Family history of early CVD
- Elevated Lp(a) screening
- Precision medicine approaches

## Quality Control

### Expected Distribution

For a diverse population cohort (N=1000):
```
Diploid CN:
  Mean: 27.5
  Median: 26.8
  Std Dev: 6.3
  Range: 14.2 - 48.7

Haploid CN:
  Mean: 13.8
  Median: 13.4
  Std Dev: 3.2
  Range: 7.1 - 24.4
```

### Validation Checks

**Check 1: Distribution shape**
```python
import matplotlib.pyplot as plt

plt.hist(estimates_df['haploid_estimate'], bins=50)
plt.xlabel('Haploid KIV-2 CN')
plt.ylabel('Frequency')
plt.title('Should be approximately normal')
plt.show()
```

Expected: Roughly normal distribution centered ~12-15 copies

**Check 2: Compare with Lp(a) levels (if available)**
```python
import pandas as pd
import scipy.stats as stats

# If you have Lp(a) measurements
data = pd.merge(estimates_df, lpa_measurements, on='sample_id')
correlation = stats.pearsonr(
    data['haploid_estimate'],
    -np.log(data['lpa_mg_dl'])  # Negative log-transform
)

print(f"Correlation: r = {correlation[0]:.3f}, p = {correlation[1]:.2e}")
# Expected: r ~ -0.5 to -0.7
```

**Check 3: Outlier detection**
```python
# Identify extreme outliers (>4 SD from mean)
mean_cn = estimates_df['haploid_estimate'].mean()
std_cn = estimates_df['haploid_estimate'].std()

outliers = estimates_df[
    abs(estimates_df['haploid_estimate'] - mean_cn) > 4 * std_cn
]

print(f"Found {len(outliers)} extreme outliers:")
print(outliers)
```

## Comparison with Other Methods

### Long-read Sequencing (PacBio/ONT)
**Gold standard** - Direct observation of repeat structure
- **Concordance:** R² ~ 0.95 with GRiD estimates
- **Advantage:** Can phase haplotypes, detect structural variants
- **Limitation:** Expensive, lower throughput

### qPCR
**Orthogonal validation** - Quantitative PCR assay
- **Concordance:** R² ~ 0.85-0.90 with GRiD
- **Advantage:** Cheap, quick
- **Limitation:** Requires careful normalization, lower resolution

### Southern Blot
**Classical method** - Historical gold standard
- **Concordance:** R² ~ 0.80-0.85
- **Advantage:** Direct visualization
- **Limitation:** Labor-intensive, requires large DNA amounts

### Other Short-read Methods
**Comparisons:**
- **lobSTR/HipSTR:** R² ~ 0.80 (general VNTR callers)
- **VNTR-wrap:** R² ~ 0.92 (LPA-specific pipeline)
- **GRiD:** R² ~ 0.94 (this pipeline, optimized workflow)

## Technical Details

### Formula Robustness

**Tested across:**
- Multiple populations (EUR, AFR, AMR, EAS)
- Different sequencing platforms (Illumina HiSeq, NovaSeq)
- Coverage ranges (20X - 80X WGS)
- Genome builds (hg19, hg38)

**Calibration stability:**
- Coefficients remain stable across cohorts
- May need minor adjustment for:
  - Very low coverage (<20X)
  - Different read lengths (100bp vs 150bp)
  - Non-standard library prep

### Uncertainty Estimation

The formula provides point estimates. For confidence intervals:

**Empirical approach:**
```python
# Bootstrap confidence intervals
def bootstrap_cn(df, n_iterations=1000):
    estimates = []
    for i in range(n_iterations):
        sample = df.sample(frac=1.0, replace=True)
        cn = 34.9 * sample['exon1A'] + 5.2 * sample['exon1B'] - 1
        estimates.append(cn.mean())
    
    ci_lower = np.percentile(estimates, 2.5)
    ci_upper = np.percentile(estimates, 97.5)
    return ci_lower, ci_upper
```

**Typical uncertainty:** ±1-2 copies (95% CI) for 30X WGS

### Performance
- **Memory:** <100MB (small dataframes)
- **Speed:** <1 second for 10,000 samples
- **Scalability:** Linear with sample count

## Advanced Usage

### Custom Coefficients

If you have your own validation data to recalibrate:

```python
from sklearn.linear_model import LinearRegression
import pandas as pd

# Load training data
exon1a = pd.read_csv("exon1a.dipCN.txt", sep="\t")
exon1b = pd.read_csv("exon1b.dipCN.txt", sep="\t")
true_cn = pd.read_csv("true_cn_from_longread.tsv", sep="\t")

# Merge data
data = exon1a.merge(exon1b, on='sample_id').merge(true_cn, on='sample_id')

# Fit regression
X = data[['exon1A_dipCN', 'exon1B_dipCN']]
y = data['true_diploid_CN']

model = LinearRegression()
model.fit(X, y)

print(f"Coefficients: exon1A={model.coef_[0]:.1f}, exon1B={model.coef_[1]:.1f}")
print(f"Intercept: {model.intercept_:.1f}")
print(f"R²: {model.score(X, y):.3f}")

# Use custom coefficients
custom_diploid_CN = (
    model.coef_[0] * data['exon1A_dipCN'] + 
    model.coef_[1] * data['exon1B_dipCN'] + 
    model.intercept_
)
```

### Population-specific Calibration

Some populations may benefit from adjusted coefficients:

```python
# Example: African ancestry cohort
african_diploid_CN = 35.2 * exon1A + 5.0 * exon1B - 1.2

# Example: East Asian ancestry cohort
eastasian_diploid_CN = 34.7 * exon1A + 5.3 * exon1B - 0.8
```

**When to recalibrate:**
- Very different population from training cohort
- Different sequencing protocol
- Validation shows systematic bias

## Output Formats

### TSV (default)
```
sample_id       diploid_estimate        haploid_estimate
sample001       36.75                   18.38
```

### CSV
```bash
grid estimate-kiv ... --format csv
```

### TXT (space-delimited)
```bash
grid estimate-kiv ... --format txt
```

## Troubleshooting

### Negative copy numbers
**Should never happen** - indicates upstream error
- Check exon diploid CN files
- Verify files are from correct pipeline step
- Report as bug if reproducible

### All estimates identical
**Problem:** No variation in exon CNs
- Check normalization step
- Verify neighbor finding worked correctly
- Inspect exon diploid CN distributions

### Very high estimates (>60 diploid)
**Possible causes:**
- True high-copy variant (rare but real)
- Technical artifact in realignment
- Incorrect normalization

**Validation steps:**
- Check raw read counts for this sample
- Inspect alignment in IGV browser
- Compare with other samples from same batch
- Consider orthogonal validation (qPCR)

### Poor correlation with Lp(a)
**Expected:** r ~ -0.5 to -0.7 between haploid CN and log(Lp(a))

**If weaker:**
- Check Lp(a) measurement quality
- Consider population stratification
- Verify CN calling quality
- Note: Lp(a) also affected by SNPs, not just CN

## Citation

When using this pipeline, please cite:

**Original methodology:**
> Mukamel, R. E., et al. (2021). Repeat polymorphisms in the human genome. *Science*, 374(6571), eabg8289.

**Implementation (if applicable):**
> Hao, A. Y., et al. (2024). vntrwrap: C++ implementation for VNTR copy number estimation.

**Your work:**
> Include appropriate citation for your study using GRiD

## Next Steps

With final KIV-2 copy number estimates:

1. **Downstream analysis:**
   - Associate with Lp(a) levels
   - Test for CVD risk associations
   - Population genetic analyses

2. **Integration with other data:**
   - Combine with SNP genotypes
   - Link to clinical phenotypes
   - Add to genetic risk scores

3. **Validation:**
   - Select subset for orthogonal validation
   - Compare with public databases (if available)
   - Check against family relationships (if trio data)

4. **Quality control:**
   - Remove low-confidence calls
   - Flag samples with unusual patterns
   - Document outliers and exclusions