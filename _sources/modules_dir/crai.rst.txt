CRAM Index Management
=====================

Overview
--------

Ensures CRAM index (.crai) files exist for efficient random access to CRAM files.

CRAM files require an index for efficient random access to specific genomic regions. This utility automatically checks for and creates the necessary `.crai` index file if it doesn't exist.

CLI Reference
-------------

.. code-block:: bash

   grid crai [OPTIONS]

**Options:**

``-c, --cram PATH``
    Input CRAM file [required]

``-r, --reference PATH``
    Reference genome FASTA [required]

``-h, --help``
    Show help message and exit

Usage Examples
--------------

Via CLI
~~~~~~~

**Basic usage:**

.. code-block:: bash

   grid crai \
       --cram sample.cram \
       --reference refs/hs37d5.fa

**Process multiple files:**

.. code-block:: bash

   # Loop through all CRAMs
   for cram in data/CRAMs/*.cram; do
       grid crai --cram "$cram" --reference refs/hs37d5.fa
   done

**In a pipeline:**

.. code-block:: bash

   # Ensure index exists before processing
   grid crai -c sample.cram -r refs/hs37d5.fa
   grid count-reads -C data/CRAMs -o counts.tsv ...

Via Python
~~~~~~~~~~

.. code-block:: python

   from grid.utils.ensure_crai import ensure_crai

   # Single file
   crai_path = ensure_crai(
       cram_path="sample.cram",
       reference="refs/hs37d5.fa"
   )
   print(f"Index ready: {crai_path}")

   # Multiple files
   from pathlib import Path
   
   cram_dir = Path("data/CRAMs")
   ref_fasta = "refs/hs37d5.fa"
   
   for cram_file in cram_dir.glob("*.cram"):
       crai_path = ensure_crai(cram_path=str(cram_file), reference=ref_fasta)
       print(f"✓ {cram_file.name}")

Description
-----------

**What it does:**

1. Checks if `.crai` index exists alongside the CRAM file
2. If missing, creates index using `samtools index`
3. Returns path to the index file
4. Validates index is properly created

**Why you need it:**

CRAM files are compressed BAM files that require an index for:

- **Random access** - Jump to specific genomic regions without reading entire file
- **Parallel processing** - Process different regions simultaneously
- **Memory efficiency** - Read only necessary data
- **Pipeline requirements** - Most downstream tools require indexed CRAMs

**Index location:**

The index is created in the same directory as the CRAM with a `.crai` extension:

.. code-block:: text

   data/
   ├── sample001.cram
   ├── sample001.cram.crai  ← Created automatically
   ├── sample002.cram
   └── sample002.cram.crai  ← Created automatically

Technical Details
-----------------

**Algorithm:**

.. code-block:: python

   def ensure_crai(cram_path, reference):
       crai_path = f"{cram_path}.crai"
       
       if not exists(crai_path):
           # Create index using samtools
           run(["samtools", "index", "-@", "4", cram_path])
       
       return crai_path

**Performance:**

- **Speed:** ~1-5 minutes per CRAM (depends on file size)
- **Memory:** <2GB RAM
- **Disk:** Index file is ~0.1% of CRAM size (e.g., 10MB for 10GB CRAM)
- **Threads:** Can specify number of threads with samtools `-@` flag

**Index structure:**

The `.crai` file contains:

- Genomic position index for rapid seeking
- Compression block offsets
- Metadata for efficient access
- Binary format (not human-readable)

Python API Documentation
------------------------

grid.utils.ensure_crai
~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: grid.utils.ensure_crai
   :members:
   :undoc-members:
   :show-inheritance:

Common Use Cases
----------------

Use Case 1: Pipeline Initialization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Ensure all CRAMs are indexed before starting analysis:

.. code-block:: bash

   #!/bin/bash
   # prepare_crams.sh
   
   CRAM_DIR="data/CRAMs"
   REF="refs/hs37d5.fa"
   
   echo "Indexing CRAMs..."
   for cram in "$CRAM_DIR"/*.cram; do
       if [ ! -f "${cram}.crai" ]; then
           echo "Indexing $(basename $cram)..."
           grid crai -c "$cram" -r "$REF"
       else
           echo "✓ $(basename $cram) already indexed"
       fi
   done
   echo "All CRAMs indexed!"

Use Case 2: Parallel Indexing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Index multiple CRAMs in parallel on HPC:

.. code-block:: bash

   #!/bin/bash
   #SBATCH --array=1-100
   #SBATCH --cpus-per-task=4
   
   # Get CRAM file for this array job
   CRAM=$(ls data/CRAMs/*.cram | sed -n ${SLURM_ARRAY_TASK_ID}p)
   
   grid crai -c "$CRAM" -r refs/hs37d5.fa

Use Case 3: Conditional Indexing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Only create index if it doesn't exist or is outdated:

.. code-block:: python

   from pathlib import Path
   from grid.utils.ensure_crai import ensure_crai
   
   def index_if_needed(cram_path, reference, force=False):
       """Index CRAM only if needed."""
       cram = Path(cram_path)
       crai = Path(f"{cram_path}.crai")
       
       # Check if index exists and is newer than CRAM
       if crai.exists() and not force:
           if crai.stat().st_mtime > cram.stat().st_mtime:
               print(f"✓ {cram.name}: Index up to date")
               return str(crai)
       
       # Create or recreate index
       print(f"Indexing {cram.name}...")
       return ensure_crai(cram_path, reference)

Dependencies
------------

**Required:**

- **samtools** - Must be installed and in PATH

  .. code-block:: bash

     # Install via conda
     conda install -c bioconda samtools
     
     # Verify installation
     samtools --version

- **pysam** - Python library for BAM/CRAM handling

  .. code-block:: bash

     pip install pysam

**Optional:**

- **parallel** - GNU parallel for batch indexing
- **slurm** - For HPC array job indexing

Troubleshooting
---------------

**Error: samtools: command not found**

.. code-block:: text

   Solution: Install samtools
   
   conda install -c bioconda samtools
   # or
   apt-get install samtools  # Ubuntu/Debian

**Error: [E::hts_open_format] fail to open file**

.. code-block:: text

   Cause: Reference genome mismatch or corrupted CRAM
   
   Solution: Verify reference matches CRAM:
   - Check CRAM header: samtools view -H sample.cram | grep @SQ
   - Ensure reference build matches (hg19 vs hg38)
   - Try redownloading CRAM if corrupted

**Error: Permission denied**

.. code-block:: text

   Cause: No write permission in CRAM directory
   
   Solution: 
   - Check directory permissions: ls -la data/CRAMs
   - Ensure you have write access
   - Or create index in writable location (not recommended)

**Warning: Index exists but empty**

.. code-block:: text

   Cause: Indexing interrupted or failed
   
   Solution: Remove and recreate
   
   rm sample.cram.crai
   grid crai -c sample.cram -r refs/hs37d5.fa

**Performance: Indexing very slow**

.. code-block:: text

   Cause: Large CRAM files or slow disk I/O
   
   Solutions:
   - Use SSD instead of HDD if possible
   - Increase samtools threads: modify ensure_crai.py
   - Process in parallel across multiple files
   - Consider pre-indexing large datasets

Quality Checks
--------------

Verify index integrity:

.. code-block:: bash

   # Check index exists and has content
   ls -lh sample.cram.crai
   
   # Verify CRAM can be accessed with index
   samtools view -H sample.cram chr6:160000000-160100000
   
   # Compare indexed vs non-indexed access time
   time samtools view sample.cram chr6:160000000-160100000 > /dev/null

Expected output:

.. code-block:: text

   -rw-r--r-- 1 user group 15M Oct 30 12:00 sample.cram.crai
   
   # Should complete in <1 second with index
   real    0m0.523s

Best Practices
--------------

1. **Index immediately** - Create indices right after downloading/generating CRAMs
2. **Store together** - Keep `.crai` files alongside `.cram` files
3. **Version control** - Recreate indices if CRAM is modified
4. **Backup strategy** - Indices are small, but can be regenerated if lost
5. **Parallel processing** - Index multiple CRAMs simultaneously to save time
6. **Check integrity** - Verify index works before starting long pipelines

Related Commands
----------------

After ensuring indices exist, proceed with analysis:

.. code-block:: bash

   # 1. Index CRAMs
   grid crai -c sample.cram -r refs/hs37d5.fa
   
   # 2. Count reads
   grid count-reads -C data/CRAMs -o counts.tsv ...
   
   # 3. Run coverage analysis
   grid mosdepth -C data/CRAMs -o coverage.tsv ...

See Also
--------

- :doc:`subset` - Subset CRAMs to specific regions
- :doc:`count-reads` - Count reads in genomic regions
- :doc:`mosdepth` - Compute coverage statistics