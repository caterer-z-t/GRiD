Read Counting
=============

**Pipeline step:** count_reads (Step 2)

**Module:** :mod:`grid.utils.count_reads`

Step 2 counts the number of reads that originate within the VNTR region of
interest for each sample. The resulting per-sample counts serve as a
locus-specific proxy for sequencing depth and feed directly into the diploid
copy number calculation in Step 4 (see :doc:`dipcn`).

----

Read Filter
-----------

For each sample, GRiD opens the CRAM/BAM file with pysam and fetches all
reads whose alignment overlaps the region :math:`[s,\, e)` on chromosome
:math:`c`. A read :math:`r` is counted if and only if it passes every one of
the following conditions simultaneously:

.. math::

   \text{flag}(r) \in \mathcal{F}
   \quad \wedge \quad
   \text{MAPQ}(r) \geq q_{\min}
   \quad \wedge \quad
   \text{chr_mate}(r) = c
   \quad \wedge \quad
   \neg\,\text{dup}(r)
   \quad \wedge \quad
   \neg\,\text{sec}(r)
   \quad \wedge \quad
   s \leq \text{pos}(r) < e

where:

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Symbol
     - Meaning
   * - :math:`\mathcal{F}`
     - User-specified set of SAM bitwise flags (``count_reads.flags`` in config).
       Each flag must be an **exact** match — not a bitmask test.
   * - :math:`q_{\min}`
     - Minimum mapping quality (``count_reads.min_mapq``).
   * - :math:`\text{chr_mate}(r) = c`
     - Both reads in the pair align to the same chromosome (proper pair).
   * - :math:`\neg\,\text{dup}(r)`
     - Read is not marked as a PCR/optical duplicate.
   * - :math:`\neg\,\text{sec}(r)`
     - Read is not a secondary alignment.
   * - :math:`s \leq \text{pos}(r) < e`
     - The read's **start position** falls within the target region boundaries.
       Reads that overlap the region but start outside it are excluded.

The per-sample read count is:

.. math::

   n_i = \bigl|\bigl\{\, r \in \mathcal{R}_i \;\mid\; \text{all conditions above}\,\bigr\}\bigr|

where :math:`\mathcal{R}_i` is the set of reads fetched for sample :math:`i`.

----

Flag Selection
--------------

The ``flags`` list controls which read orientations are counted. The filter
uses an **exact** flag comparison (not a bitwise AND), so each value in the
list targets a specific orientation. The recommended defaults capture
properly paired reads on both strands:

.. list-table::
   :header-rows: 1
   :widths: 10 90

   * - Flag
     - Orientation
   * - ``83``
     - Proper pair, read on reverse strand, mate on forward strand (R1 reverse)
   * - ``147``
     - Proper pair, mate on reverse strand, read on forward strand (R2 reverse)
   * - ``81``
     - Read on reverse strand, mate on forward strand (not requiring proper pair flag)
   * - ``145``
     - Mate on reverse strand (not requiring proper pair flag)

Flags can be added or removed to match the read orientation profile of a
specific library preparation. Use the
`Broad Picard flag explainer <https://broadinstitute.github.io/picard/explain-flags.html>`_
to decode or construct flag values.

----

Implementation
--------------

Samples are processed in parallel using a ``ThreadPoolExecutor`` with the
thread count set by ``threads`` in the config. Results are written to the
output file incrementally as each sample completes, using a threading lock
for safe concurrent writes.

The output is a tab-separated file with one row per sample:

.. code-block:: text

   HG00096.hg38.cram    1842
   HG00097.hg38.cram    1956
   NA12878.hg38.cram    1703
   ...

----

Downstream Use
--------------

The count :math:`n_i` is used in Step 4 (:doc:`dipcn`) as the numerator of
the diploid CN estimator. For each sample :math:`i`, its locus read count is
scaled relative to the mean count of its nearest neighbors :math:`\mathcal{N}_i`:

.. math::

   \widehat{\text{dipCN}}_i
   = 2 \cdot \frac{n_i / \bar{d}_i}{\frac{1}{|\mathcal{N}_i|}\sum_{j \in \mathcal{N}_i} n_j / \bar{d}_j}

where :math:`\bar{d}_i` is the genome-wide depth scale for sample :math:`i`
derived from the mosdepth normalization step. See :doc:`dipcn` for the full
derivation.
