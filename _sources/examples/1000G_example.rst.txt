1000 Genomes Full Pipeline Example
==================================

This example demonstrates Stage 3 of the GRiD workflow: dynamically generating ``config.yaml`` and running the complete GRiD pipeline on the 1000 Genomes Project dataset for the *LPA* KIV-2 locus (GRCh38).

The full workflow comprises three sequential modular scripts:

1. :doc:`prep_data_example`: Downloads inputs, streams CRAM slices, and writes locus coordinates.
2. :doc:`IBS_example`: *(Optional)* Computes IBS neighbors for haploid copy number estimation.
3. **GRiD Execution:** Generates configuration and runs GRiD steps (indexing, read counting, depth normalization, and genotype calling).

Stage 3 Execution Script
------------------------

``03-run-grid.sh`` checks for the presence of ``data/ibs_neighbors_chr6.tsv.gz`` from Stage 2. If present, haploid estimation (``compute_haploid_genotypes.run``) is set to ``True``. Otherwise, it falls back to diploid-only mode.

.. literalinclude:: ../../../examples/03-run-grid.sh
   :language: bash
   :linenos:
   :caption: 03-run-grid.sh — Stage 3: Config generation and pipeline execution

Full Pipeline Usage
-------------------

Run all three stages sequentially from your working directory:

.. code-block:: bash

   # Step 1: Prepare data and stream CRAMs
   sbatch examples/01-prep-data.sh

   # Step 2: (Optional) Compute IBS neighbors for haploid CN
   sbatch examples/02-run-ibs.sh

   # Step 3: Generate config and run GRiD
   sbatch examples/03-run-grid.sh