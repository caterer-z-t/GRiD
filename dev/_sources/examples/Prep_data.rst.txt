Data Preparation Example
========================

This example covers Stage 1 of the GRiD workflow: fetching reference inputs, parsing locus coordinates, streaming regional CRAM slices from 1000 Genomes, and building sample manifests.

It sets up the necessary data directory hierarchy and exports coordinate variables to ``data/pipeline_vars.sh`` so downstream steps (IBS computation and GRiD configuration) operate on identical locus bounds.

.. literalinclude:: ../../../examples/01-prep-data.sh
   :language: bash
   :linenos:
   :caption: 01-prep-data.sh — Stage 1: Reference, CRAM streaming, and manifest prep

What This Script Does
---------------------

1. **Locus Coordinates:** Downloads the coding VNTR regions file, extracts *LPA* coordinates, and writes ``data/pipeline_vars.sh``.
2. **Reference Genome:** Downloads and indexes the GRCh38 no-alt reference genome via ``samtools faidx``.
3. **Sample Manifest:** Fetches the 1000 Genomes panel file and extracts sample IDs.
4. **CRAM Streaming:** Parallel-streams *LPA* regional slices directly from EBI using ``samtools view`` and indexes them locally into ``crams/``.
5. **GRiD Manifest:** Writes ``data/grid_samples.txt`` containing valid, non-empty sample IDs for processing.

Usage
-----

Submit to SLURM:

.. code-block:: bash

   sbatch examples/01-prep-data.sh

Or execute locally:

.. code-block:: bash

   bash examples/01-prep-data.sh