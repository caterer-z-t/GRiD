Coverage Normalization
======================

**Pipeline step:** normalize_mosdepth (Step 4)

**Module:** :mod:`grid.utils.normalize_mosdepth`

Raw depth values vary between samples due to differences in sequencing depth,
library preparation, and batch effects. Normalization removes these technical
sources of variation so that differences in the depth matrix reflect true
copy number variation rather than measurement artefacts.

Depth Matrix
------------

Let :math:`N` be the number of samples and :math:`M` the number of genomic bins
after filtering. Define the raw depth matrix:

.. math::

   D \in \mathbb{R}^{N \times M}, \qquad
   D_{ij} =
   \begin{cases}
   d_{ij} & \text{if sample } i \text{ has depth for bin } j \\
   \mathrm{NaN} & \text{otherwise}
   \end{cases}

Step 1 — Within-Individual Normalization
-----------------------------------------

Each sample's depth values are divided by that sample's mean depth across all
non-missing bins, removing global sequencing depth differences:

.. math::

   \bar{D}_i =
   \frac{1}{|\{j : D_{ij} \neq \mathrm{NaN}\}|}
   \sum_{\substack{j=1 \\ D_{ij} \neq \mathrm{NaN}}}^{M} D_{ij}

.. math::

   D^{(1)}_{ij} =
   \begin{cases}
   D_{ij} \,/\, \bar{D}_i & \bar{D}_i > 0 \\
   \mathrm{NaN} & \bar{D}_i = 0
   \end{cases}

Step 2 — Across-Individual Normalization
-----------------------------------------

For each bin :math:`j`, compute the population mean and variance of the
within-normalized depths:

.. math::

   \mu_j =
   \frac{1}{|\{i : D^{(1)}_{ij} \neq \mathrm{NaN}\}|}
   \sum_{\substack{i=1 \\ D^{(1)}_{ij} \neq \mathrm{NaN}}}^{N} D^{(1)}_{ij}

.. math::

   \sigma_j^{2} =
   \frac{1}{|\{i : D^{(1)}_{ij} \neq \mathrm{NaN}\}|}
   \sum_{\substack{i=1 \\ D^{(1)}_{ij} \neq \mathrm{NaN}}}^{N}
   \left(D^{(1)}_{ij} - \mu_j\right)^{2}

The final normalized matrix uses a Poisson-inspired z-score transformation
(dividing by :math:`\sqrt{\mu_j}` rather than :math:`\sigma_j`) to stabilize
variance across bins with different mean depths:

.. math::

   D^{\mathrm{norm}}_{ij} =
   \begin{cases}
   \dfrac{D^{(1)}_{ij} - \mu_j}{\sqrt{\mu_j}} & \mu_j > 0 \\[0.8em]
   \mathrm{NaN} & \mu_j \leq 0
   \end{cases}

Step 3 — Variance Ratio Filtering
-----------------------------------

To focus downstream neighbor search on the most informative regions
(those with true copy number variation), bins are ranked by their
variance ratio:

.. math::

   \mathrm{VarRatio}_j =
   \begin{cases}
   \sigma_j^{2} \,/\, \mu_j & \mu_j > 0 \\
   \mathrm{NaN} & \mu_j \leq 0
   \end{cases}

Only the top fraction (default: top 10%) of bins by variance ratio are
retained for the nearest-neighbor step. This preferentially selects bins
in and around the VNTR region where copy number variation is greatest.

**Output:** A gzip-compressed matrix of shape :math:`N \times M_{\mathrm{filtered}}`
containing :math:`D^{\mathrm{norm}}` values for the retained bins.
