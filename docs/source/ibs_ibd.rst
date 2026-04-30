Computing IBS/IBD Neighbors
=======================================================

Step 7 of GRiD (haplotype inference) requires a file of identity-by-descent (IBD) /
identity-by-state (IBS) haplotype neighbors as input. This file is produced by
``computeIBSpbwt.cpp``, a C++ program from Hujoel et al. (2026) that identifies
the closest haplotype matches for each individual using the Positional Burrows-Wheeler
Transform (PBWT).

This page covers how to obtain the source, compile it, prepare the required inputs,
and run it to generate the neighbors file used by ``HI_inference``.

.. note::

   ``computeIBSpbwt.cpp`` is not included in the GRiD repository. It is part of the
   supplementary code released with Hujoel et al. (2026) Nature
   (`doi:10.1038/s41586-025-09886-z <https://doi.org/10.1038/s41586-025-09886-z>`_).
   The source is available directly at:
   `mhujoel/STRs — computeIBSpbwt.cpp <https://github.com/mhujoel/STRs/blob/main/cpp_files/computeIBSpbwt.cpp>`_.

----

Dependencies
------------

``computeIBSpbwt.cpp`` depends on:

1. **Eagle v2.4.1 header files** — specifically ``Types.hpp``, ``FileUtils.cpp``,
   ``StringUtils.cpp``, ``MemoryUtils.cpp``, ``Timer.cpp``, and ``MapInterpolater.cpp``.
   Download Eagle from:

   .. code-block:: bash

      wget https://alkesgroup.broadinstitute.org/Eagle/downloads/Eagle_v2.4.1.tar.gz
      tar -xzf Eagle_v2.4.1.tar.gz
      # headers are in Eagle_v2.4.1/src/

2. **Boost iostreams** (≥ 1.58.0):

   .. code-block:: bash

      # via conda
      conda install -c conda-forge boost

      # or download manually
      wget https://archives.boost.org/release/1.58.0/source/boost_1_58_0.tar.gz

3. **zlib** — usually available on HPC systems; install if missing:

   .. code-block:: bash

      conda install -c conda-forge zlib

4. **OpenMP** — included with GCC; no separate install needed.

----

Compilation
-----------

Replace ``/path/to/Eagle/src``, ``/path/to/boost/include``, and
``/path/to/boost/lib`` with your actual paths.

**Standard build:**

.. code-block:: bash

   g++ -O2 -fopenmp -Wall computeIBSpbwt.cpp -o computeIBSpbwt \
       -I/path/to/Eagle_v2.4.1/src \
       -I/path/to/boost/include \
       -L/path/to/boost/lib \
       -Wl,-Bstatic -lboost_iostreams -Wl,-Bdynamic -lz

**Static build (recommended for HPC clusters):**

.. code-block:: bash

   g++ -O2 -fopenmp -Wall -static-libgcc -static-libstdc++ computeIBSpbwt.cpp -o computeIBSpbwt \
       -I/path/to/Eagle_v2.4.1/src \
       -I/path/to/boost/include \
       -L/path/to/boost/lib \
       -Wl,-Bstatic -lboost_iostreams -lz

Verify the build:

.. code-block:: bash

   ./computeIBSpbwt
   # should print: Usage: - arg1 = chr (integer) ...

----

Required Input Files
--------------------

The program takes 8 arguments in order:

.. list-table::
   :header-rows: 1
   :widths: 10 20 70

   * - Arg
     - Name
     - Description
   * - 1
     - ``chr``
     - Chromosome number as integer (e.g. ``6`` for *LPA*)
   * - 2
     - ``focal_bp``
     - Focal base pair position in the same genome build as the BGEN file. For *LPA* KIV-2 in hg19, use the center of the VNTR region (:math:`\sim` ``160690000``).
   * - 3
     - ``bgen_file``
     - Phased genotype data in BGEN v1.2 format (layout 2, compressed, phased). See below.
   * - 4
     - ``sample_file``
     - Oxford-format sample file corresponding to the BGEN file. See below.
   * - 5
     - ``genetic_map_file``
     - Genetic map in cM units, in the same genome build as the BGEN file. See below.
   * - 6
     - ``num_neighbors``
     - Number of IBS neighbors to find per haplotype. Recommended: ``200``
   * - 7
     - ``threads``
     - Number of CPU threads for parallel computation.
   * - 8
     - ``output_file``
     - Output file path. Use a ``.gz`` extension for automatic gzip compression.

Phased BGEN File
~~~~~~~~~~~~~~~~

The BGEN file must contain **phased** SNP-array or WGS genotypes for the chromosome
of interest. Requirements:

- BGEN v1.2, layout 2, compressed (``CompressedSNPBlocks=1``)
- Phased (``Phased=1``)
- 16-bit probability encoding (``bgenBits=16``)
- Biallelic variants only (``K=2``)

Sample File
~~~~~~~~~~~

Oxford-format ``.sample`` file with two header lines followed by one row per individual:

.. code-block:: text

   ID_1 ID_2 missing
   0 0 0
   SAMPLE001 SAMPLE001 0
   SAMPLE002 SAMPLE002 0
   ...

The program uses column ``ID_1`` as the sample identifier. The identifier must match
the IDs in the diploid CN file passed to ``hi_inference``.

Genetic Map File
~~~~~~~~~~~~~~~~

A tab-separated file with columns: ``position``, ``COMBINED_rate``, ``Genetic_Map(cM)``.
Download pre-built maps from Eagle's website:

.. code-block:: bash

   wget https://alkesgroup.broadinstitute.org/Eagle/downloads/tables/genetic_map_hg19.tar.gz
   # or for hg38:
   wget https://alkesgroup.broadinstitute.org/Eagle/downloads/tables/genetic_map_hg38.tar.gz

Running the Program
-------------------

Example command for *LPA* KIV-2 (hg19, chr6):

.. code-block:: bash

   ./computeIBSpbwt \
       6 \
       160690000 \
       phased_chr6.bgen \
       phased_chr6.sample \
       genetic_map_chr6_combined_b37.txt \
       200 \
       16 \
       ibs_neighbors_chr6.tsv.gz

On a SLURM cluster:

.. code-block:: bash

   #!/bin/bash
   #SBATCH --job-name=computeIBS
   #SBATCH --cpus-per-task=16
   #SBATCH --mem=64G
   #SBATCH --time=4:00:00

   ./computeIBSpbwt 6 160690000 \
       phased_chr6.bgen phased_chr6.sample \
       genetic_map_chr6_combined_b37.txt \
       200 16 \
       ibs_neighbors_chr6.tsv.gz

Output Format
-------------

The output is a tab-separated file with one row per haplotype-neighbor pair:

.. code-block:: text

   ID      hap  nbrInd  cMlen   cMedge  IDnbr      hapNbr
   10001   1    1       8.45    4.22    20034      2
   10001   1    2       7.13    3.56    30017      1
   10001   2    1       9.01    5.10    20034      1
   ...

.. list-table::
   :header-rows: 1
   :widths: 15 85

   * - Column
     - Description
   * - ``ID``
     - Sample identifier
   * - ``hap``
     - Haplotype index (1 or 2)
   * - ``nbrInd``
     - Neighbor rank (1 = closest)
   * - ``cMlen``
     - Total IBS match length in cM (weighted combination of left and right extents)
   * - ``cMedge``
     - Minimum of left and right IBS match edge lengths in cM
   * - ``IDnbr``
     - Neighbor sample identifier
   * - ``hapNbr``
     - Neighbor haplotype index (1 or 2)

This file is passed directly to GRiD Step 7 via ``compute_haploid_genotypes.ibs_output``
in the config.

----

Passing the Output to GRiD
---------------------------

Once the neighbors file is generated, point to it in your config:

.. code-block:: yaml

   compute_haploid_genotypes:
     run: True
     output_file_prefix: "haploid_genotypes"
     phased_vcf: "path/to/phased.vcf.gz"
     ibs_output: "path/to/ibs_neighbors_chr6.tsv.gz"
     min_neighbors: 1
     max_neighbors: 10
     n_iters: 100

See :doc:`algorithms/hi_inference` for details on how the neighbors file is used
in the iterative phasing algorithm.
