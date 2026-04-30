Haplotype Inference
===================

Decomposes diploid copy number estimates into haplotype-specific values using an
iterative IBD-based phasing algorithm. For each sample, haplotype copy numbers
are updated proportionally to the mean of IBD-matched haplotype neighbors until
convergence. This step is locus-agnostic — it requires only diploid CN estimates
and a phased IBD neighbors file.

.. automodule:: grid.utils.hi_inference
   :members:
   :undoc-members:
   :show-inheritance:
