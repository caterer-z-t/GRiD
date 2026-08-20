IBS/IBD Neighbor Computation Example
====================================

This example demonstrates Stage 2 of the GRiD workflow: using ``computeIBSpbwt`` to generate the IBS neighbors file required by GRiD for **haploid copy number estimation** (Step 5).

The script automatically sources locus coordinates (chromosome and focal midpoint) derived during Stage 1 from ``data/pipeline_vars.sh``.

.. literalinclude:: ../../../examples/02-run-ibs.sh
   :language: bash
   :linenos:
   :caption: 02-run-ibs.sh — Stage 2: Compute IBS neighbors via PBWT

Key Parameters & Environment Variables
--------------------------------------

The coordinates and variables used by ``02-run-ibs.sh`` are sourced dynamically from Stage 1 (``data/pipeline_vars.sh``) or can be passed via command line options:

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Parameter / Variable
     - Description
   * - ``--num-neighbors``
     - Number of IBS neighbors per haplotype (default: ``200``).
   * - ``CHR``
     - Target chromosome (e.g., ``6``, exported by Stage 1).
   * - ``FOCAL_BP``
     - Base pair position at center of VNTR (e.g., *LPA* KIV-2 midpoint in hg38, exported by Stage 1).
   * - ``COMPUTE_IBS``
     - Path to the auto-downloaded or local ``computeIBSpbwt`` binary.
   * - ``PHASED_VCF_URL``
     - 1000 Genomes high-coverage phased VCF URL streamed for the target region via ``tabix``.
   * - ``GENETIC_MAP``
     - Eagle hg38 genetic map auto-downloaded and filtered for chromosome ``CHR``.
   * - ``OUTPUT_FILE``
     - Compressed output path (``data/ibs_neighbors_chr${CHR}.tsv.gz``).

Usage
-----

Submit to SLURM after ``01-prep-data.sh`` completes:

.. code-block:: bash

   sbatch examples/02-run-ibs.sh

Override the default number of neighbors:

.. code-block:: bash

   sbatch examples/02-run-ibs.sh --num-neighbors 100

Passing Output to GRiD
----------------------

Stage 3 (``03-run-grid.sh``) automatically detects ``data/ibs_neighbors_chr${CHR}.tsv.gz`` upon completion and embeds the path directly into ``config.yaml``:

.. code-block:: yaml

   compute_haploid_genotypes:
     run: True
     output_file_prefix: "haploid_genotypes"
     ibs_output: "data/ibs_neighbors_chr6.tsv.gz"
     min_neighbors: 1
     max_neighbors: 10
     n_iters: 100
     method: "ibs"