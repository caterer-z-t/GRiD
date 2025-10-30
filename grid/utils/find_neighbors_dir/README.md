# Neighbor Finding Module

Identify nearest neighbors among samples based on genome-wide normalized coverage patterns. Used for sample-specific normalization in copy number estimation.

## Purpose

Copy number estimation benefits from comparing each sample to genetically similar individuals rather than the entire cohort. This module finds the most similar samples based on coverage profiles, which:

1. **Reduces batch effects** - Similar samples likely sequenced together
2. **Improves accuracy** - Better baseline for CN estimation  
3. **Handles population structure** - Accounts for ancestry-related variation
4. **Minimizes noise** - Related samples have similar technical characteristics

## Main Module

### `find_neighbors.py`
Located in parent directory: `grid/utils/find_neighbors.py`

**Key Function:**
```python
def find_neighbors(
    input_file: str,
    output_file: str,
    zmax: float = 2.0,
    n_neighbors: int = 500,
    sigma2_max: float = 1000.0
) -> None
```

**Parameters:**
- `input_file`: Normalized z-depth file from normalize_mosdepth (`.gz`)
- `output_file`: Output neighbor file (auto-suffixed with `.zMax{zmax}.txt.gz`)
- `zmax`: Maximum z-score clipping threshold (default: 2.0)
- `n_neighbors`: Number of nearest neighbors to find (default: 500)
- `sigma2_max`: Maximum variance threshold for region filtering (default: 1000.0)

**Output Format:**
```
sample_id       neighbor_1      neighbor_2      neighbor_3      ...
sample001       sample042       sample103       sample205       ...
sample002       sample073       sample192       sample301       ...
```

Each row contains a sample and its N nearest neighbors, ranked by similarity.

## Usage

### Via CLI
```bash
grid find-neighbors \
    --input-file output/normalized.zdepth.gz \
    --output-file output/neighbors.txt \
    --zmax 2.0 \
    --n-neighbors 500 \
    --sigma2-max 1000.0
```

The output will be automatically named: `output/neighbors.zMax2.0.txt.gz`

### As Python Module
```python
from grid.utils.find_neighbors import find_neighbors

find_neighbors(
    input_file="output/normalized.zdepth.gz",
    output_file="output/neighbors.txt",
    zmax=2.0,
    n_neighbors=500,
    sigma2_max=1000.0
)
```

## Algorithm Details

### Step 1: Load and Filter Data
1. **Load z-depths:** Read normalized coverage matrix
2. **Clip outliers:** Cap z-scores at ±`zmax` to reduce outlier influence
3. **Filter high-variance regions:** Remove positions with variance > `sigma2_max`

**Z-clipping example:**
```
Original:  [-3.2, -0.5, 1.2, 5.8, 0.3]
After clip (zmax=2): [-2.0, -0.5, 1.2, 2.0, 0.3]
```

### Step 2: Compute Pairwise Distances
Calculate Euclidean distance between all sample pairs:

```
distance(A, B) = sqrt(Σ(z_A[i] - z_B[i])²)
```

Where the sum is over all retained genomic positions.

**Similarity interpretation:**
- **Small distance** = Similar coverage patterns = Good neighbors
- **Large distance** = Different patterns = Poor neighbors

### Step 3: Rank Neighbors
For each sample:
1. Sort all other samples by distance (ascending)
2. Select top N closest samples as neighbors
3. Write neighbor list to output file

## Parameter Tuning

### `zmax` (default: 2.0)
**Purpose:** Reduce influence of extreme outliers

- **Lower (1.5):** More conservative, less outlier influence
- **Higher (3.0):** Keep more extreme variation
- **Recommended:** 2.0 for most applications

**Effect:**
- Too low: May lose real biological signal
- Too high: Outliers dominate distance calculation

### `n_neighbors` (default: 500)
**Purpose:** Number of similar samples to identify

- **Fewer (100-200):** Very stringent, only closest matches
- **More (500-1000):** Include more distant but still similar samples
- **Recommended:**
  - Small cohort (<500): 100-200 neighbors
  - Medium cohort (500-2000): 200-500 neighbors  
  - Large cohort (>2000): 500-1000 neighbors

**Trade-off:**
- Too few: May miss good reference samples
- Too many: Include dissimilar samples, reduce specificity

### `sigma2_max` (default: 1000.0)
**Purpose:** Filter extremely variable regions

**Effect:**
- Lower (500): More stringent filtering, fewer positions
- Higher (2000): Keep more variable positions
- Infinite: Keep all positions

**Recommended:**
- Standard VNTR: 1000
- Highly variable CNVs: 2000
- Conserved regions: 500

## Output Interpretation

### Neighbor Ranking
Neighbors are ordered by similarity:
```
sample001: neighbor1 (closest), neighbor2, neighbor3, ... neighbor500 (farthest)
```

**Using neighbor data:**
The top neighbors are used in `compute_dipcn` for sample-specific normalization:
```python
# Typically use top 100-200 neighbors for CN estimation
top_neighbors = neighbor_list[:200]
reference_coverage = mean(top_neighbors.coverage)
```

### Quality Metrics

**Good neighborhood:**
- Neighbors show consistent coverage patterns in VNTR
- Similar sequencing depth across genome
- Likely from same sequencing batch/population

**Poor neighborhood:**
- High variance among neighbors
- May indicate population outlier
- Could suggest unique structural variants

## Technical Details

### Computational Complexity
- **Memory:** O(samples × positions) - full matrix in memory
- **Time:** O(samples² × positions) - pairwise comparisons
- **Parallelization:** Not currently implemented (TODO)

**Performance estimates:**
- 1,000 samples, 1,000 positions: ~1-2 minutes
- 5,000 samples, 1,000 positions: ~10-20 minutes
- 10,000 samples: Memory intensive (>32GB RAM)

### Memory Optimization
For very large cohorts (>5000 samples):
1. Use chunked processing
2. Implement approximate nearest neighbor search (e.g., FAISS)
3. Reduce number of positions used

## Use in Pipeline

The neighbor file is critical for downstream analysis:

```
normalize_mosdepth → find_neighbors → compute_dipcn → estimate_kiv
                            ↓
                   Neighbor-based CN normalization
```

**In compute_dipcn:**
- Load neighbor lists
- For each sample, extract coverage of top N neighbors
- Use neighbor median/mean as reference baseline
- Compute sample CN relative to neighbor baseline

## Validation

### Manual Check
Inspect a few samples' neighbors:
```python
import pandas as pd

neighbors = pd.read_csv("neighbors.zMax2.0.txt.gz", sep="\t")
print(neighbors.iloc[0])  # Check first sample's neighbors

# Are neighbors from same population/batch?
# Do they have similar coverage characteristics?
```

### PCA Visualization
```python
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt

# Project samples into 2D
pca = PCA(n_components=2)
coords = pca.fit_transform(zdepth_matrix)

plt.scatter(coords[:, 0], coords[:, 1])
plt.xlabel("PC1")
plt.ylabel("PC2")
plt.title("Sample similarity space")
```

Neighbors should cluster together in PCA space.

## Troubleshooting

**All distances similar:**
- May indicate insufficient variation in data
- Check that normalization worked correctly
- Verify region contains informative positions

**Extreme distances:**
- Check for outlier samples with very unusual coverage
- May indicate technical failures or true large CNVs
- Consider removing outliers before neighbor finding

**Memory errors:**
- Reduce `n_neighbors`
- Process samples in batches
- Use sparse matrix representation (TODO)