1000 Genomes Full Pipeline Example
=====================================

This example runs the complete GRiD pipeline (from raw CRAMs through diploid
copy number estimation) using the publicly available 1000 Genomes Project
high-coverage WGS dataset. No login or data access agreement is required.

.. literalinclude:: ../../../examples/1000G_example.sh
   :language: bash
   :linenos:

Usage
-----

Run all samples from a single superpopulation:

.. code-block:: bash

   bash examples/1000G_example.sh --pop EUR

Limit to 50 samples for a quick test:

.. code-block:: bash

   bash examples/1000G_example.sh --pop EUR --n 50

Submit to SLURM:

.. code-block:: bash

   sbatch --cpus-per-task=16 --mem=32G examples/1000G_example.sh --pop EUR --n 100

Override the working directory:

.. code-block:: bash

   WORK_DIR=/scratch/my_run bash examples/1000G_example.sh --pop AFR

Next Step — Haplotype Inference
---------------------------------

The pipeline above estimates **diploid** copy number. To decompose into
haplotype-specific copy numbers, you need an IBS neighbors file first.
See :doc:`IBS_example` for instructions, then re-run the pipeline with
``compute_haploid_genotypes.run: True`` in the generated config.
