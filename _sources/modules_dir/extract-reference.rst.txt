Extract Reference Sequences
===========================

Extract FASTA Sequences from Reference
---------------------------------------

Extract specific genomic regions from a reference genome based on BED coordinates.

Python API Documentation
------------------------

grid.utils.extract_reference
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: grid.utils.extract_reference
   :members:
   :undoc-members:
   :show-inheritance:

Usage
-----

Via CLI
~~~~~~~

.. code-block:: bash

   grid extract-reference \
       --reference-fa refs/hs37d5.fa \
       --bed-file lpa_regions.bed \
       --output-dir refs/ \
       --output-prefix lpa_kiv2

Via Python
~~~~~~~~~~

.. code-block:: python

   from grid.utils.extract_reference import extract_reference

   extract_reference(
       reference_fa="refs/hs37d5.fa",
       bed_file="lpa_regions.bed",
       output_dir="refs/",
       output_prefix="lpa_kiv2"
   )

Description
-----------

This utility extracts sequences from a reference genome for specific regions:

1. Reads BED file with target coordinates
2. Extracts sequences using samtools faidx
3. Writes to new FASTA file
4. Creates index (.fai) for new reference

Use cases:

- **LPA Reference** - Extract KIV-2 region for realignment
- **Custom References** - Build VNTR-specific references
- **Validation** - Create test references for development

BED File Format
---------------

Standard BED format (0-based coordinates):

.. code-block:: text

   chr6    160500000    160510000    KIV2_region1
   chr6    160520000    160530000    KIV2_region2

Columns:

1. Chromosome name
2. Start position (0-based)
3. End position (exclusive)
4. Region name (optional)

Output
------

Creates:

- ``<prefix>.fa`` - Extracted sequences in FASTA format
- ``<prefix>.fa.fai`` - Index file for the extracted reference

The output FASTA will have headers based on the BED regions:

.. code-block:: text

   >chr6:160500000-160510000
   ATCGATCGATCG...
   >chr6:160520000-160530000
   GCTAGCTAGCTA...

Dependencies
------------

- **samtools** - Must be installed and in PATH (uses faidx)
- **pysam** - Python library for FASTA handling

Notes
-----

- Reference genome must be indexed (``.fai`` file required)
- Chromosome names in BED must match reference exactly
- Extracted regions maintain original coordinates in headers