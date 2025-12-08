# Mosdepth
1. check for crai file
2. build mosdepth command
3. run mosdepth
4. wait for mosdepth output
5. compute coverage
    $$ 
    \text{Coverage} =
    \begin{cases}
    \displaystyle 
    \left\lfloor 
    100 \cdot 
    \frac{
    \sum\limits_{i \in R} \left( \overline{C}_i \cdot 
    \bigl( \min(e, r_{i,e}) - \max(s, r_{i,s}) \bigr) \right)
    }{
    \sum\limits_{i \in R} 
    \bigl( \min(e, r_{i,e}) - \max(s, r_{i,s}) \bigr)
    }
    \right\rceil,
    & \text{if } \sum\limits_{i \in R} \text{overlap}_i > 0, \\[1.3em]
    0, & \text{otherwise}.
    \end{cases}
    $$

# Normalize Mosdepth
1. check / collect individuals (extension `.mosdepth.global.dist.txt`)
2. load repeat mask
3. Process one individual
    - a. extract regions 
    - b. if region is none, then find surrounding 1000 bp region
4. build depth matrix from regions
    - Combine all individuals’ regions into a unified set of genomic intervals:
    $$
    \mathcal{R} = \{ (s_1, e_1), (s_2, e_2), \ldots, (s_M, e_M) \}.
    $$

    - Construct a matrix
    $$
    D \in \mathbb{R}^{N \times M},
    $$

    Where,
    $$
    D_{ij} =
    \begin{cases}
    d_{ij}, & \text{if sample } i \text{ has depth for region } j, \\[0.5em]
    \mathrm{NaN}, & \text{otherwise}.
    \end{cases}
    $$
    ### Input: Depth values per individual

    | Sample | Region (start, end) | Depth |
    |--------|----------------------|-------|
    | A      | (100, 200)           | 12.5  |
    | A      | (300, 400)           | 18.0  |
    | B      | (100, 200)           | 10.1  |
    | B      | (500, 600)           | 22.0  |
    | C      | (300, 400)           | 17.5  |

    ### Depth Matrix (rows = samples, columns = regions)

    | Sample | (100,200) | (300,400) | (500,600) |
    |--------|-----------|-----------|-----------|
    | A      | 12.5      | 18.0      | NaN       |
    | B      | 10.1      | NaN       | 22.0      |
    | C      | NaN       | 17.5      | NaN       |

5. Normalize matrix
    - a. within-individual
    $$
    \bar{D}_i
    =
    \frac{1}{\left| \{ j : D_{ij} \neq \mathrm{NaN} \} \right|}
    \sum_{\substack{j=1 \\ D_{ij} \neq \mathrm{NaN}}}^{M}
    D_{ij}, 
    \\
    \text{Where}, 
    \\
    D^{(1)}_{ij}
    =
    \begin{cases}
    \dfrac{D_{ij}}{\bar{D}_i}, & \bar{D}_i > 0, \\[0.8em]
    \mathrm{NaN}, & \bar{D}_i = 0.
    \end{cases}
    $$
    - b. across-individual
    $$
    \mu_j
    =
    \frac{1}{\left| \{ i : D^{(1)}_{ij} \neq \mathrm{NaN} \} \right|}
    \sum_{\substack{i=1 \\ D^{(1)}_{ij} \neq \mathrm{NaN}}}^{N} 
    D^{(1)}_{ij},
    \\
    \sigma_j^{2}
    =
    \frac{1}{\left| \{ i : D^{(1)}_{ij} \neq \mathrm{NaN} \} \right|}
    \sum_{\substack{i=1 \\ D^{(1)}_{ij} \neq \mathrm{NaN}}}^{N}
    \left(D^{(1)}_{ij} - \mu_j\right)^{2}
    $$
    - c. variance ratio
    $$ 
    \mathrm{VarRatio}_j =
    \begin{cases}
    \dfrac{\sigma_j^{2}}{\mu_j}, & \mu_j > 0, \\[0.8em]
    \mathrm{NaN}, & \mu_j \le 0.
    \end{cases}
    $$
    - d. final matrix
    $$
    D^{\mathrm{norm}}_{ij}
    =
    \begin{cases}
    \dfrac{D^{(1)}_{ij} - \mu_j}{\sqrt{\mu_j}}, & \mu_j > 0, \\[0.8em]
    \mathrm{NaN}, & \mu_j \le 0.
    \end{cases}
    $$
6. select high variance regions
    - select top 10% of variance regions (aka looking for vntr region)

# Find neighbors
1. read normalized depth matrix 
2. filter regions by variance
    - a. calculate variance ratio for each region
    - b. find regions that pass variance and fraction ratio
3. find neighbors sklearn
    - a. Input data
    $$
    D \in \mathcal{R}^{M \times N}, \\
    - M = \text{number of regions} \\
    - N = \text{number of individuals} \\
    - D_{ij} = \text{normalized depth for region } i \text{ and } j \\
    \text{transpose the matrix to} \\
    X = D^T \in \mathcal{R}^{N \times M}, X_{k\ell} = D_{k\ell}
    $$

    - b. Clip and replace NaNs
    $$
    X^`_{k\ell} = \begin{cases}
    z_{\text{max}}, &  X_{k\ell} >z_{\text{max}}\\
    -z_{\text{max}}, & X_{k\ell} <-z_{\text{max}}\\
    0, & X_{k\ell} = \text{NaN} \\
    X_{k\ell}, & \text{otherwise}
    \end{cases}
    $$

    - c. calculate eucledian distance between individuals
    $$
    d_{ij} = ||X_i^` - X_j^`||_2 = \sqrt{\sum_{\ell=1}^M (X_{i\ell}^` - X_{j\ell}^`)^2}
    $$

    - d. find nearst neighbors (using scikit-learn)
    $$
    N_i = \Bigl\{\, j \in \{1, \dots, N\} \setminus \{i\} \;\Big|\; d_{ij} \text{ is among the smallest } n_{\text{neighbors}} \text{ distances} \Bigr\}, \\
    \text{where } d_{ij} = \| X_i - X_j \|_2 \text{ and } i = 1, \dots, N.
    $$
4. save neighbors

# Extract Reference Command
1. extract specific region from fasta + bed file
2. save this reference region as a fasta file

# 