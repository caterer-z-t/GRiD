Haplotype Inference
===================

**Pipeline step:** compute_haploid_genotypes (Step 5)

**Module:** :mod:`grid.utils.hi_inference`

Diploid copy number estimates from Step 4 (denoted IRR — individual repeat ratio)
cannot be decomposed into per-haplotype contributions without phase information.
GRiD infers haplotype-specific copy numbers iteratively by identifying
haplotype-matched neighbors, following Hujoel et al. (2026).

Two neighbor sources are supported via the ``method`` config parameter:

- **IBS** (``method: "ibs"``) — neighbors from ``computeIBSpbwt`` PBWT matching output.
  Each neighbor is weighted equally.
- **IBD** (``method: "ibd"``) — neighbors from iLASH IBD segment output.
  Set ``weighted: True`` to apply a Lorentzian distance+match weight to each segment.

For more information on the neighbor file formats, see :doc:`../ibs_ibd`.

.. note::

   This step is **locus-agnostic** — it operates on diploid copy number estimates
   and a haplotype neighbors file, with no assumptions about which locus produced
   those estimates. It was developed and validated on *LPA* KIV-2, but can be applied
   to any VNTR locus for which phased neighbor data are available.

Notation
--------

For :math:`N` individuals, define:

- :math:`\mathrm{IRR}_i` — diploid copy number estimate for individual :math:`i`
- :math:`h_{i,k}` — haplotype copy number for individual :math:`i`, haplotype :math:`k \in \{1, 2\}`
- :math:`\mathcal{H}_{i,k}` — set of IBD-matched haplotype neighbors for haplotype :math:`k` of individual :math:`i`

The constraint :math:`h_{i,1} + h_{i,2} = \mathrm{IRR}_i` must hold at every iteration.

Initialization
--------------

Each haplotype is initialized to half the diploid estimate:

.. math::

   h_{i,1}^{(0)} = h_{i,2}^{(0)} = \frac{\mathrm{IRR}_i}{2}

Individuals with fewer than :math:`\mathrm{MIN\_NBR}` IBD neighbors on either
haplotype are excluded from iterative phasing and receive mean-imputed values.

Iterative Update
----------------

At each iteration :math:`t`, a weighted mean haplotype copy number across
matched neighbors is computed for each haplotype:

.. math::

   \bar{h}^{(t)}_{\mathcal{H}_{i,k}} =
   \frac{\displaystyle\sum_{(j,\, k') \in \mathcal{H}_{i,k}} w_{j,k'} \cdot h^{(t)}_{j,\, k'}}
        {\displaystyle\sum_{(j,\, k') \in \mathcal{H}_{i,k}} w_{j,k'}}

where :math:`(j, k')` denotes haplotype :math:`k'` of neighbor :math:`j` matched
to haplotype :math:`k` of individual :math:`i`.

**IBS method:** all weights :math:`w_{j,k'} = 1`.

**IBD method (unweighted):** all weights :math:`w_{j,k'} = 1`.

**IBD method (weighted):** the weight accounts for segment proximity to the
target region and IBD match quality:

.. math::

   w_{j,k'} = \frac{\lambda}{\delta_{j,k'} + \lambda} \cdot m_{j,k'}

where :math:`\delta_{j,k'}` is the physical distance (bp) from the IBD segment
to the target region (0 if overlapping), :math:`m_{j,k'}` is the iLASH match
score, and :math:`\lambda` is the ``weight_scale`` parameter (default 1,000,000 bp).

Each haplotype is updated proportionally, preserving the diploid constraint:

.. math::

   h^{(t+1)}_{i,k} =
   \mathrm{IRR}_i \cdot
   \frac{\bar{h}^{(t)}_{\mathcal{H}_{i,k}}}
        {\bar{h}^{(t)}_{\mathcal{H}_{i,1}} + \bar{h}^{(t)}_{\mathcal{H}_{i,2}}}

This redistributes the diploid total between haplotypes in proportion to their
neighbors' copy numbers. In practice the algorithm runs for a fixed number of
iterations (default: :math:`T = 100`).

Mean Imputation for Unphased Samples
--------------------------------------

For individuals with insufficient IBD neighbors, haplotype copy numbers
are mean-imputed using the cohort mean on each haplotype:

.. math::

   h^{\mathrm{imp}}_{i,k} = \bar{h}^{(T)}_{\mathcal{H}_{i,k}}

where the mean is taken over any available neighbors (or the cohort mean
if no neighbors exist).

Output
------

The output file contains one row per individual with six columns:

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Column
     - Description
   * - ``ID``
     - Sample identifier
   * - ``IRRs``
     - Original diploid copy number estimate
   * - ``hap1phased``
     - Inferred haplotype 1 copy number (iterative phasing)
   * - ``hap2phased``
     - Inferred haplotype 2 copy number (iterative phasing)
   * - ``hap1imp``
     - Mean-imputed haplotype 1 copy number (neighbor mean)
   * - ``hap2imp``
     - Mean-imputed haplotype 2 copy number (neighbor mean)
