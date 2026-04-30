Nearest-Neighbor Identification
================================

**Pipeline step:** find_neighbors (Step 5)

**Module:** :mod:`grid.utils.find_neighbors`

Copy number normalization requires a reference set of individuals with similar
genome-wide sequencing characteristics. Rather than using a fixed external panel,
GRiD finds each individual's nearest neighbors within the cohort itself using
the normalized depth matrix. This naturally controls for ancestry, batch effects,
and sequencing platform differences without requiring explicit covariates.

Input
-----

Let :math:`D^{\mathrm{norm}} \in \mathbb{R}^{N \times M}` be the filtered normalized
depth matrix from the previous step, where :math:`N` is the number of individuals
and :math:`M` is the number of retained high-variance bins. Define the transposed
input matrix:

.. math::

   X = \left(D^{\mathrm{norm}}\right)^{\top} \in \mathbb{R}^{N \times M},
   \qquad X_{k\ell} = D^{\mathrm{norm}}_{k\ell}

so that each row :math:`X_k` represents one individual's depth profile.

Step 1 — Z-score Clipping and NaN Replacement
-----------------------------------------------

Extreme depth values at outlier bins are clipped to prevent them from dominating
distance calculations. Missing values are set to zero (the population mean in
normalized space):

.. math::

   X'_{k\ell} =
   \begin{cases}
   z_{\max}  & X_{k\ell} > z_{\max} \\
   -z_{\max} & X_{k\ell} < -z_{\max} \\
   0         & X_{k\ell} = \mathrm{NaN} \\
   X_{k\ell} & \text{otherwise}
   \end{cases}

where :math:`z_{\max}` is configurable (default: 2.0).

Step 2 — Euclidean Distance
-----------------------------

Pairwise distances between individuals are computed in the clipped normalized
depth space:

.. math::

   d_{ij} = \left\| X'_i - X'_j \right\|_2
           = \sqrt{\sum_{\ell=1}^{M} \left(X'_{i\ell} - X'_{j\ell}\right)^2}

Step 3 — k-Nearest Neighbor Selection
---------------------------------------

For each individual :math:`i`, the :math:`k` nearest neighbors are those with
the smallest Euclidean distances (excluding the individual itself):

.. math::

   \mathcal{N}_i = \Bigl\{
     j \in \{1, \ldots, N\} \setminus \{i\}
     \;\Big|\;
     d_{ij} \text{ is among the } k \text{ smallest}
   \Bigr\}

This is computed efficiently using scikit-learn's
``NearestNeighbors`` with a ball-tree index. The default is :math:`k = 500`
neighbors, of which the top 200 are used in the diploid CN step.

**Output:** A gzip-compressed file listing, for each individual, their neighbor
IDs, neighbor depth scales, and normalized distances. The depth scale per
individual (their mean mosdepth coverage) is also stored for use in the
diploid CN computation.
