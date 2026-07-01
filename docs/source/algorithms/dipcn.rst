Diploid Copy Number Estimation
===============================

**Pipeline step:** compute_dipcn (Step 4)

**Module:** :mod:`grid.utils.compute_dipcn`

Raw read counts in the target VNTR region depend on both the true copy number and
each sample's sequencing depth. To isolate the copy number signal, GRiD
normalizes each sample's counts relative to its nearest neighbors — individuals
with similar genome-wide depth profiles (and therefore similar technical
characteristics). This step is fully locus-agnostic.

Notation
--------

For sample :math:`i`:

- :math:`r_i` — read count in the target VNTR region
- :math:`s_i` — depth scale (mean mosdepth coverage, stored in the neighbors file)
- :math:`\mathcal{N}_i` — set of top-:math:`k` nearest neighbors (default :math:`k = 200`)

For each neighbor :math:`j \in \mathcal{N}_i`:

- :math:`r_j` — neighbor's read count
- :math:`s_j` — neighbor's depth scale

Depth-Scaled Normalization
---------------------------

Each read count is first divided by its sample's depth scale to remove the effect
of sequencing depth:

.. math::

   \tilde{r}_i = \frac{r_i}{s_i}

The neighbor mean of depth-scaled counts provides a cohort-matched reference:

.. math::

   \bar{r}_{\mathcal{N}_i} = \frac{1}{|\mathcal{N}_i|}
   \sum_{j \in \mathcal{N}_i} \frac{r_j}{s_j}

Diploid Copy Number Estimate
-----------------------------

The diploid copy number estimate is the ratio of the sample's depth-scaled count
to the neighbor mean:

.. math::

   \widehat{\mathrm{dipCN}}_i =
   \frac{\tilde{r}_i}{\bar{r}_{\mathcal{N}_i}} =
   \frac{r_i / s_i}{\dfrac{1}{|\mathcal{N}_i|} \displaystyle\sum_{j \in \mathcal{N}_i} r_j / s_j}

Under the assumption that the average neighbor has diploid copy number equal to
the cohort mean, :math:`\widehat{\mathrm{dipCN}}_i` is a relative estimate of
copy number normalized to the cohort. When used with the LPA realignment step,
this is computed separately for each exon classification (e.g. 1A, 1B-KIV2, 1B-KIV3);
for other loci a single count per sample is sufficient.

**Output:** A TSV file per exon type with sample ID and normalized diploid copy
number estimate.
