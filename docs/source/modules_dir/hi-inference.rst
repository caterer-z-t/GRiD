Haplotype Inference
===================

Decomposes diploid copy number estimates into haplotype-specific values using
neighbor-based phasing. Supports two neighbor sources (IBS via ``computeIBSpbwt``,
or IBD via iLASH) and five phasing algorithms selectable via config.
This step is locus-agnostic — it requires only diploid CN estimates and a phased
neighbors file.

For a full description of the algorithms, neighbor quality-control options,
and configuration reference, see :doc:`../algorithms/hi_inference`.

.. automodule:: grid.utils.hi_inference
   :members:
   :undoc-members:
   :show-inheritance: