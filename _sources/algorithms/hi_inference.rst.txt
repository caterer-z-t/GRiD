Haplotype Inference
===================

**Pipeline step:** compute_haploid_genotypes (Step 5)

**Module:** :mod:`grid.utils.hi_inference`

Diploid copy number estimates from Step 4 (denoted IRR — individual repeat ratio)
cannot be decomposed into per-haplotype contributions without phase information.
GRiD infers haplotype-specific copy numbers by identifying haplotype-matched
neighbors and iteratively (or directly) partitioning each individual's diploid
total between their two haplotypes.

Two neighbor sources are supported via the ``method`` config parameter:

- **IBS** (``method: "ibs"``) — neighbors from ``computeIBSpbwt`` PBWT matching output.
- **IBD** (``method: "ibd"``) — neighbors from iLASH IBD segment output.

Five phasing algorithms are available via the ``phasing_method`` config parameter,
described in detail below.

For neighbor file formats and how to generate them, see :doc:`../ibs_ibd`.

.. note::

   This step is **locus-agnostic** — it operates on diploid copy number estimates
   and a haplotype neighbors file, with no assumptions about which locus produced
   those estimates. It was developed and validated on *LPA* KIV-2, but can be applied
   to any VNTR locus for which phased neighbor data are available.

----

Notation
--------

For :math:`N` individuals, define:

- :math:`\mathrm{IRR}_i` — diploid copy number estimate for individual :math:`i`
- :math:`h_{i,k}^{(t)}` — haplotype copy number for individual :math:`i`,
  haplotype :math:`k \in \{0, 1\}`, at iteration :math:`t`
- :math:`\mathcal{H}_{i,k}` — set of haplotype neighbors for haplotype :math:`k`
  of individual :math:`i`, each entry :math:`(j, k', w)` giving neighbor individual
  :math:`j`, neighbor haplotype :math:`k'`, and weight :math:`w`
- :math:`\bar{h}^{(t)}_{i,k}` — aggregated neighbor signal for haplotype :math:`k`
  of individual :math:`i` at iteration :math:`t`

The diploid constraint :math:`h_{i,0} + h_{i,1} = \mathrm{IRR}_i` is preserved
at every iteration for all phasing methods.

----

Initialization
--------------

All methods initialize each haplotype to half the diploid estimate:

.. math::

   h_{i,0}^{(0)} = h_{i,1}^{(0)} = \frac{\mathrm{IRR}_i}{2}

Only individuals with at least :math:`\mathrm{MIN\_NBR}` neighbors on **both**
haplotypes are phased; others receive mean-imputed values (see `Mean Imputation for Unphased Samples`_).

----

Neighbor Aggregation Functions
-------------------------------

All iterative methods compute a neighbor aggregate :math:`\bar{h}^{(t)}_{i,k}`
from the current haplotype values of neighbors in :math:`\mathcal{H}_{i,k}`.
The aggregation function is controlled by the ``agg_fn`` config parameter:

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - ``agg_fn``
     - Description
   * - ``mean`` *(default)*
     - Arithmetic mean of neighbor haplotype values, ignoring weights.
   * - ``median``
     - Median of neighbor haplotype values; robust to copy-number outliers.
   * - ``weighted_cM``
     - Weighted mean using segment length (cM) as weight. For IBS neighbors,
       all weights are 1 and this reduces to ``mean``.
   * - ``weighted_score``
     - Weighted mean using iLASH match score as weight.
   * - ``weighted_cM_score``
     - Weighted mean using the product of segment length and match score as weight.
   * - ``trimmed_mean``
     - 10 % trimmed mean (trims the top and bottom decile before averaging).

For IBS neighbors, all segment weights are 1 so ``mean``, ``weighted_cM``,
``weighted_score``, and ``weighted_cM_score`` are all equivalent.

----

Phasing Algorithms
------------------

.. _phasing-iterative:

``iterative`` — Original Hujoel Algorithm
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The default algorithm, ported directly from Hujoel et al. (2026).
At each iteration the diploid total is partitioned in proportion to the
ratio of each haplotype's neighbor aggregate:

.. math::

   \bar{h}^{(t)}_{i,k} = \mathrm{agg\_fn}\!\left(\left\{
       h^{(t)}_{j,\, k'} \;\middle|\; (j, k', w) \in \mathcal{H}_{i,k},\;
       h^{(t)}_{j,k'} \text{ not NaN}
   \right\}\right)

.. math::

   h^{(t+1)}_{i,k} =
   \mathrm{IRR}_i \cdot
   \frac{\bar{h}^{(t)}_{i,k}}
        {\bar{h}^{(t)}_{i,0} + \bar{h}^{(t)}_{i,1}}

**Config:** ``phasing_method: iterative`` |
**Tunes:** ``agg_fn``, ``n_iters``

.. _phasing-graph-prop:

``graph_prop`` — Graph Propagation *(recommended)*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A label-spreading variant that blends each individual's own prior signal
(diploid/2) with the neighbor aggregate at every iteration, controlled by
a mixing parameter :math:`\alpha`:

.. math::

   \tilde{h}^{(t)}_{i,k} =
   \alpha \cdot \frac{\mathrm{IRR}_i}{2}
   \;+\; (1 - \alpha) \cdot \bar{h}^{(t)}_{i,k}

.. math::

   h^{(t+1)}_{i,k} =
   \mathrm{IRR}_i \cdot
   \frac{\tilde{h}^{(t)}_{i,k}}
        {\tilde{h}^{(t)}_{i,0} + \tilde{h}^{(t)}_{i,1}}

The :math:`\alpha` parameter controls how strongly each individual's own
diploid signal anchors the estimate. Higher :math:`\alpha` (e.g. 0.7)
trusts the individual's own signal more and is more resistant to being
pulled by noisy neighbors. Lower :math:`\alpha` relies more heavily on
the neighborhood.

.. note::

   Empirically, ``graph_prop`` with ``alpha: 0.7``, ``min_match: 0.95``,
   ``min_length: 2.0``, and ``max_neighbors: 20`` achieved the best
   correlation with Bionano optical mapping ground truth on *LPA* KIV-2
   (Pearson *r* = 0.83 short allele, 0.92 long allele).

**Config:** ``phasing_method: graph_prop`` |
**Tunes:** ``alpha``, ``n_iters``

``graph_prop_weighted``
~~~~~~~~~~~~~~~~~~~~~~~~

Identical to ``graph_prop`` but uses the configured ``agg_fn`` (weighted
mean, median, etc.) for the neighbor aggregate instead of plain mean.
Useful when segment-length or match-score weighting is expected to improve
discrimination.

**Config:** ``phasing_method: graph_prop_weighted`` |
**Tunes:** ``alpha``, ``agg_fn``, ``n_iters``

``oneshot`` — Single-Pass
~~~~~~~~~~~~~~~~~~~~~~~~~~

A non-iterative baseline. Neighbor values are fixed at the seed
(:math:`\mathrm{IRR}_i / 2`) and each sample's estimate is computed
independently in a single pass. Fast and useful as a lower bound on
iterative performance.

.. math::

   h_{i,k}^{\mathrm{oneshot}} =
   \mathrm{IRR}_i \cdot
   \frac{\mathrm{agg\_fn}\!\left(\left\{
       \tfrac{\mathrm{IRR}_j}{2} \;\middle|\; (j, k', w) \in \mathcal{H}_{i,k}
   \right\}\right)}
   {\sum_{k=0}^{1}
       \mathrm{agg\_fn}\!\left(\left\{
           \tfrac{\mathrm{IRR}_j}{2} \;\middle|\; (j, k', w) \in \mathcal{H}_{i,k}
       \right\}\right)}

**Config:** ``phasing_method: oneshot`` |
**Tunes:** ``agg_fn``

``rank_based`` — Rank-Normalised Iterative
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Like ``iterative``, but before each aggregation step all haplotype values
are rank-normalised to :math:`[0, 1]` across the cohort. This reduces the
influence of extreme copy-number outliers on the neighborhood signal.

.. math::

   r^{(t)}_{i,k} = \frac{\mathrm{rank}\!\left(h^{(t)}_{i,k}\right)}{|\text{phased haplotypes}|}

Aggregation and partitioning proceed on the rank-normalised values; the
diploid total is then used to recover the original scale.

**Config:** ``phasing_method: rank_based`` |
**Tunes:** ``agg_fn``, ``n_iters``

----

IBD Neighbor Loading and Quality Control
-----------------------------------------

When ``method: ibd``, neighbors are loaded from an iLASH ``.match`` file.
Two additional quality-control mechanisms beyond a flat match-score threshold
are available.

.. _tiered-neighbors:

Two-Tier Neighbor Priority
~~~~~~~~~~~~~~~~~~~~~~~~~~

When ``start_bp`` and ``end_bp`` are set in the top-level config (the VNTR
coordinates), IBD segments are assigned to one of two tiers before filling
each haplotype's neighbor list:

**Tier 1 (gold):** segment physically overlaps the target VNTR region
:math:`[\mathrm{start\_bp},\, \mathrm{end\_bp}]` **and** passes the quality gate.
These neighbors share IBD directly at the locus of interest and are the most
informative for haplotype-specific CN estimation.

**Tier 2 (fill):** segment passes the quality gate but does not overlap the VNTR
region. Used only to fill any remaining slots in a haplotype's neighbor list after
Tier 1 neighbors are exhausted.

Within each tier, segments are ranked by :math:`\mathrm{cM} \times \mathrm{score}`
descending. Set ``use_tiers: False`` to disable tiering and treat all qualifying
segments equally.

.. _adaptive-threshold:

Length-Adaptive Match Threshold
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Short IBD segments have a higher false-positive rate than long segments
and should require a stricter match score. Setting ``match_at_min_length``
enables a length-dependent threshold that scales continuously from a
strict ceiling at ``min_length`` cM down to the lenient floor ``min_match``
at long segments:

.. math::

   \theta(L) = \mathrm{min\_match}
   \;+\; (\mathrm{match\_at\_min\_length} - \mathrm{min\_match})
   \cdot \frac{\mathrm{min\_length}}{L}

At :math:`L = \mathrm{min\_length}`, :math:`\theta = \mathrm{match\_at\_min\_length}`.
As :math:`L \to \infty`, :math:`\theta \to \mathrm{min\_match}`.
If ``match_at_min_length`` is not set, a flat threshold of ``min_match`` is
applied to all segments.

----

Mean Imputation for Unphased Samples
--------------------------------------

Individuals with fewer than ``min_neighbors`` neighbors on either haplotype
cannot be phased iteratively. For these individuals, haplotype copy numbers
are mean-imputed from whatever neighbors are available:

.. math::

   h^{\mathrm{imp}}_{i,k} =
   \frac{\displaystyle\sum_{(j,k',w) \in \mathcal{H}_{i,k}} w \cdot h^{(T)}_{j,k'}}
        {\displaystyle\sum_{(j,k',w) \in \mathcal{H}_{i,k}} w}

If no neighbors exist, the cohort mean :math:`\bar{\mathrm{IRR}} / 2` is used.
Imputed values are written to ``hap1imp`` / ``hap2imp`` in the output file.

----

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
     - Inferred haplotype 0 copy number (iterative phasing); ``NaN`` if sample
       had insufficient neighbors
   * - ``hap2phased``
     - Inferred haplotype 1 copy number (iterative phasing); ``NaN`` if sample
       had insufficient neighbors
   * - ``hap1imp``
     - Mean-imputed haplotype 0 copy number (neighbor mean or cohort mean)
   * - ``hap2imp``
     - Mean-imputed haplotype 1 copy number (neighbor mean or cohort mean)

----

Configuration Reference
------------------------

The full set of ``compute_haploid_genotypes`` config keys:

.. list-table::
   :header-rows: 1
   :widths: 30 12 58

   * - Key
     - Default
     - Description
   * - ``method``
     - ``"ibs"``
     - Neighbor source: ``"ibs"`` or ``"ibd"``.
   * - ``phasing_method``
     - ``"graph_prop"``
     - Phasing algorithm: ``"iterative"``, ``"graph_prop"``,
       ``"graph_prop_weighted"``, ``"oneshot"``, or ``"rank_based"``.
   * - ``agg_fn``
     - ``"mean"``
     - Neighbor aggregation: ``"mean"``, ``"median"``, ``"weighted_cM"``,
       ``"weighted_score"``, ``"weighted_cM_score"``, or ``"trimmed_mean"``.
       Ignored by plain ``graph_prop``.
   * - ``alpha``
     - ``0.7``
     - Mixing weight for ``graph_prop`` / ``graph_prop_weighted``.
       Range (0, 1); higher values trust the individual's own diploid signal more.
   * - ``n_iters``
     - ``20``
     - Number of iterations. Ignored by ``oneshot``.
   * - ``min_neighbors``
     - ``1``
     - Minimum neighbors required on both haplotypes to phase a sample.
   * - ``max_neighbors``
     - ``20``
     - Maximum neighbors loaded per haplotype.
   * - ``ibs_output``
     - —
     - Path to IBS neighbors file. Required when ``method: ibs``.
   * - ``ibd_output``
     - —
     - Path to iLASH ``.match`` file. Required when ``method: ibd``.
   * - ``min_length``
     - ``0.5``
     - Minimum IBD segment length in cM.
   * - ``min_match``
     - ``0.70``
     - Minimum iLASH match score (lenient floor for long segments when
       ``match_at_min_length`` is set).
   * - ``match_at_min_length``
     - ``None``
     - Strict match score required at ``min_length`` cM. Enables the
       length-adaptive threshold. Leave unset for a flat threshold.
   * - ``weighted``
     - ``False``
     - Apply :math:`\mathrm{match}^3 \times \mathrm{cM}` weighting to IBD
       neighbor weights during loading.
   * - ``weight_scale``
     - ``1000000``
     - :math:`\lambda` parameter (bp) for Lorentzian distance weighting
       (only when ``weighted: True``).

**Example config (IBD, recommended settings):**

.. code-block:: yaml

   compute_haploid_genotypes:
     run: True
     output_file_prefix: "haploid_genotypes"
     method: "ibd"
     phasing_method: "graph_prop"
     ibd_output: "path/to/output.match"
     alpha: 0.7
     n_iters: 20
     max_neighbors: 20
     min_length: 2.0
     min_match: 0.95

**Example config (IBS, original algorithm):**

.. code-block:: yaml

   compute_haploid_genotypes:
     run: True
     output_file_prefix: "haploid_genotypes"
     method: "ibs"
     phasing_method: "iterative"
     ibs_output: "path/to/ibs_neighbors_chr6.tsv.gz"
     min_neighbors: 1
     max_neighbors: 10
     n_iters: 20