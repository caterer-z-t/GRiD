IBS/IBD Neighbor Computation Example
=====================================

This example shows how to run ``computeIBSpbwt`` to produce the IBS neighbors file
required by GRiD Step 5 (haplotype inference). The script targets *LPA* KIV-2 on
chromosome 6 in hg19, but the same pattern applies to any locus — adjust
``CHR``, ``FOCAL_BP``, and the matching BGEN/genetic-map files accordingly.

Before running, compile ``computeIBSpbwt.cpp`` and obtain the required input files.
See the :doc:`../ibs_ibd` page for full instructions.

.. literalinclude:: ../../../examples/IBS_example.sh
   :language: bash
   :linenos:
   :caption: IBS_example.sh — SLURM/local script for computing IBS neighbors

Key Parameters
--------------

Edit the ``Paths`` and ``Parameters`` sections at the top of the script before submitting:

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Variable
     - Description
   * - ``COMPUTE_IBS``
     - Path to the compiled ``computeIBSpbwt`` binary
   * - ``BGEN_FILE``
     - Phased BGEN v1.2 file for the chromosome of interest
   * - ``SAMPLE_FILE``
     - Oxford-format ``.sample`` file matching the BGEN
   * - ``GENETIC_MAP``
     - Eagle genetic map for the chromosome (download from Eagle website)
   * - ``OUTPUT_FILE``
     - Output path; use ``.tsv.gz`` for automatic gzip compression
   * - ``FOCAL_BP``
     - Base pair position at the center of the VNTR region (hg19: ``160690000``)
   * - ``NUM_NEIGHBORS``
     - Number of IBS neighbors per haplotype (recommended: ``200``)
   * - ``THREADS``
     - CPU threads; match to ``--cpus-per-task`` in the SLURM header

Usage
-----

Submit to SLURM:

.. code-block:: bash

   sbatch IBS_example.sh

Or run locally:

.. code-block:: bash

   bash IBS_example.sh

Passing Output to GRiD
-----------------------

Once complete, point to the output file in your GRiD config:

.. code-block:: yaml

   compute_haploid_genotypes:
     run: True
     output_file_prefix: "haploid_genotypes"
     ibs_output: "path/to/ibs_neighbors_chr6.tsv.gz"
     min_neighbors: 1
     max_neighbors: 10
     n_iters: 100
