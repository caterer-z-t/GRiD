Haplotype Inference
===================

**Pipeline step:** compute_haploid_genotypes (Step 5)

**Module:** :mod:`grid.utils.hi_inference`

Diploid copy number estimates from Step 4 (denoted IRR — individual repeat ratio)
cannot be decomposed into per-haplotype contributions without phase information.
GRiD uses identity-by-descent (IBD) segments to identify haplotype-matched neighbors
and infers haplotype-specific copy numbers iteratively, following Hujoel et al. (2026).
For more information on the IBD neighbors file, see :doc:`../ibs_ibd`.

.. note::

   This step is **locus-agnostic** — it operates on diploid copy number estimates
   and a phased IBD neighbors file, with no assumptions about which locus produced
   those estimates. It was developed and validated on *LPA* KIV-2, but can be applied
   to any VNTR locus for which a phased VCF and IBD output are available.

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

At each iteration :math:`t`, the mean haplotype copy number across IBD-matched
neighbors is computed for each haplotype:

.. math::

   \bar{h}^{(t)}_{\mathcal{H}_{i,k}} =
   \frac{1}{|\mathcal{H}_{i,k}|}
   \sum_{(j,\, k') \in \mathcal{H}_{i,k}} h^{(t)}_{j,\, k'}

where :math:`(j, k')` denotes haplotype :math:`k'` of neighbor :math:`j` that
shares an IBD segment with haplotype :math:`k` of individual :math:`i`.

Each haplotype is then updated proportionally, preserving the diploid constraint:

.. math::

   h^{(t+1)}_{i,k} =
   \mathrm{IRR}_i \cdot
   \frac{\bar{h}^{(t)}_{\mathcal{H}_{i,k}}}
        {\bar{h}^{(t)}_{\mathcal{H}_{i,1}} + \bar{h}^{(t)}_{\mathcal{H}_{i,2}}}

This update redistributes the diploid total between the two haplotypes in
proportion to how much copy number their IBD neighbors carry. Convergence
is reached when haplotype values stabilize; in practice the algorithm runs
for a fixed number of iterations (default: :math:`T = 100`).

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
