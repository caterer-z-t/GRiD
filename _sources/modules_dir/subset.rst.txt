CRAM Subsetting
===============

Subset CRAM to Genomic Region
------------------------------

Extract reads from a specific genomic region to create a smaller, focused CRAM file.

Python API Documentation
------------------------

grid.utils.subset_cram
~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: grid.utils.subset_cram
   :members:
   :undoc-members:
   :show-inheritance:

Usage
-----

Via CLI
~~~~~~~

.. code-block:: bash

   # Using region string
   grid subset \
       --cram input.cram \
       --region chr6:160000000-160100000 \
       --output subset.cram \
       --reference refs/hs37d5.fa

   # Using separate coordinates
   grid subset \
       --cram input.cram \
       --chrom chr6 \
       --start 160000000 \
       --end 160100000 \
       --output subset.cram \
       --reference refs/hs37d5.fa

Via Python
~~~~~~~~~~

.. code-block:: python

   from grid.utils.subset_cram import subset_cram

   subset_path = subset_cram(
       cram_path="input.cram",
       region="chr6:160000000-160100000",
       output_path="subset.cram",
       reference="refs/hs37d5.fa"
   )

Description
-----------

This utility extracts reads overlapping a specific genomic region:

1. Uses samtools view to extract region
2. Maintains CRAM format and compression
3. Creates index (.crai) for output file
4. Preserves header information

Use cases:

- **Testing** - Work with smaller files during development
- **Sharing** - Extract specific regions for collaboration
- **Performance** - Faster processing on region of interest
- **Storage** - Reduce file sizes for targeted analysis

Region Format
-------------

Regions can be specified as:

- **Full format**: ``chr6:160000000-160100000``
- **Chromosome only**: ``chr6`` (entire chromosome)
- **No end**: ``chr6:160000000`` (from position to end)

Coordinates are 0-based, following SAM/BAM conventions.

Dependencies
------------

- **samtools** - Must be installed and in PATH
- **pysam** - Python library for BAM/CRAM handling

Output
------

Creates two files:

- ``<output>.cram`` - Subset CRAM file
- ``<output>.cram.crai`` - Index for subset CRAM