Coverage Estimation
===================

**Pipeline step:** mosdepth (Step 3)

**Module:** :mod:`grid.utils.mosdepth`

GRiD uses `mosdepth <https://github.com/brentp/mosdepth>`_ to compute read depth
across the genome in fixed-width bins (default 1,000 bp). Because the target VNTR region
boundaries rarely align exactly to bin edges, a single bin's depth is not sufficient.
Instead, GRiD computes an overlap-weighted mean depth across all bins that intersect
the target region. This step is fully locus-agnostic — the genomic coordinates
are set in the config file.

Overlap-Weighted Mean Depth
----------------------------

Let :math:`[s, e]` denote the genomic region of interest (set via ``chrom``, ``start_bp``, ``end_bp`` in the config),
and let :math:`R` be the set of mosdepth bins that overlap this region.
For each bin :math:`i \in R`, let :math:`\bar{C}_i` be its mean read depth,
:math:`r_{i,s}` its start coordinate, and :math:`r_{i,e}` its end coordinate.

The overlap length of bin :math:`i` with the target region is:

.. math::

   \ell_i = \min(e,\, r_{i,e}) - \max(s,\, r_{i,s})

The overlap-weighted mean depth is then:

.. math::

   \text{Coverage} =
   \begin{cases}
   \left\lfloor
   100 \cdot
   \dfrac{\displaystyle\sum_{i \in R} \bar{C}_i \cdot \ell_i}
         {\displaystyle\sum_{i \in R} \ell_i}
   \right\rceil
   & \text{if } \displaystyle\sum_{i \in R} \ell_i > 0 \\[1.5em]
   0 & \text{otherwise}
   \end{cases}

The result is multiplied by 100 and rounded to an integer to preserve two decimal
places of precision in integer storage.

**Output:** A TSV file with one row per sample containing the sample ID and its
integer-encoded coverage value.
