Getting Started
===============

GRiD is driven by a single YAML configuration file. This page walks through every
section — what each parameter does, why it is needed, and how to obtain or prepare
the required input files.

Copy the bundled example and edit it for your cohort:

.. code-block:: bash

   cp grid/example_config.yaml my_config.yaml
   # edit my_config.yaml, then:
   grid wgs my_config.yaml

----

Global Settings
---------------

These fields apply to every pipeline step.

.. code-block:: yaml

   samples_file: "path/to/sample_ids.txt"
   directory_loc: "path/to/crams/"
   reference_genome: "path/to/hg38.fa"
   output_dir: "path/to/output/"
   threads: 8
   file_type: "cram"
   chrom: "chr6"
   start_bp: 160605062
   end_bp: 160647661
   output_file_type: "tsv"

.. list-table::
   :header-rows: 1
   :widths: 22 12 66

   * - Parameter
     - Type
     - Description
   * - ``samples_file``
     - path
     - Plain-text file with one sample ID per line. IDs must match the file stems
       of the CRAM/BAM files in ``directory_loc``.
   * - ``directory_loc``
     - path
     - Directory containing the CRAM (or BAM) files for your cohort.
   * - ``reference_genome``
     - path
     - Reference genome FASTA (hg19 or hg38). Must be indexed (``samtools faidx``).
       Download from `UCSC <https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/>`_,
       the `Broad Institute <https://s3.amazonaws.com/broad-references/broad-references-readme.html>`_, or
       `1000 Genomes <https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/>`_.
   * - ``output_dir``
     - path
     - Directory where all pipeline outputs are written. Created if it does not exist.
   * - ``threads``
     - int
     - Number of CPU threads shared across steps. Match to your SLURM
       ``--cpus-per-task`` allocation.
   * - ``file_type``
     - str
     - Sequencing file format: ``"cram"`` or ``"bam"``.
   * - ``chrom``
     - str
     - Target chromosome in the same notation as your CRAM headers
       (e.g. ``"chr6"`` for hg38, ``"6"`` for hg19).
   * - ``start_bp``
     - int
     - Start coordinate (0-based) of the VNTR region of interest.
       For *LPA* KIV-2 in hg38: ``160605062``.
   * - ``end_bp``
     - int
     - End coordinate of the VNTR region of interest.
       For *LPA* KIV-2 in hg38: ``160647661``.
   * - ``output_file_type``
     - str
     - Output table format: ``"tsv"`` or ``"csv"``.

.. note::

   **Supported VNTR loci.** Mukamel et al. (2021) *Science*
   (`doi:10.1126/science.abg8289 <https://www.science.org/doi/10.1126/science.abg8289>`_)
   catalogued **734 coding VNTRs** with elevated IBD2 rates across the human genome,
   identifying loci where copy number is likely to vary and be heritable. Their full code
   is available at `Zenodo 4776804 <https://zenodo.org/record/4776804>`_.

   The same file is bundled with GRiD at 
   `734_possible_coding_vntr_regions.IBD2R_gt_0.25.uniq.txt <https://github.com/caterer-z-t/GRiD/tree/main/files/734_possible_coding_vntr_regions.IBD2R_gt_0.25.uniq.txt>`_. 
   To run GRiD on a different locus, simply update the ``chrom``, ``start_bp``, and ``end_bp`` fields in your 
   config file to match the coordinates of your target region. For example, to run on *LPA* KIV-2.

====

Step 1 — ``index``
-------------------

Verifies or creates ``.crai`` / ``.bai`` index files for all CRAM/BAM files in
``directory_loc``. Only needed if your files are not already indexed.

.. code-block:: yaml

   index:
     run: False
     output_file_prefix: "index_file_results"

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Parameter
     - Description
   * - ``run``
     - Set ``True`` to run this step. Skip if files are already indexed.
   * - ``output_file_prefix``
     - Prefix for the log output file.

Step 2 — ``count_reads``
-------------------------

Counts properly paired reads overlapping the target VNTR region from each CRAM/BAM
file. The resulting per-sample read counts are used downstream as a proxy for
locus-specific sequencing depth.

.. code-block:: yaml

   count_reads:
     run: True
     output_file_prefix: "read_counts"
     min_mapq: 1
     flags:
       - 83    # proper pair, read reverse strand
       - 147   # proper pair, mate reverse strand
       - 81    # read reverse strand
       - 145   # mate reverse strand

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Parameter
     - Description
   * - ``run``
     - Set ``True`` to run this step.
   * - ``output_file_prefix``
     - Prefix for the output read-count file.
   * - ``min_mapq``
     - Minimum mapping quality. Reads below this threshold are ignored.
       ``1`` excludes unmapped reads while retaining multi-mappers; higher
       values (e.g. ``20``) restrict to uniquely mapped reads.
   * - ``flags``
     - List of SAM bitwise flags to include. Only reads whose flag exactly
       matches one of the listed values are counted. The defaults capture
       properly paired reads on both strands. See
       `SAM flag documentation <https://broadinstitute.github.io/picard/explain-flags.html>`_
       for flag definitions.

Step 3 — ``mosdepth``
----------------------

Runs `mosdepth <https://github.com/brentp/mosdepth>`_ in windowed mode to compute
binned read depth across the full genome for each sample. The resulting depth
profiles are used for normalization and neighbor finding.

``mosdepth`` must be available on your ``PATH``:

.. code-block:: bash

   conda install -c bioconda mosdepth
   # or: module load mosdepth  (HPC)

.. code-block:: yaml

   mosdepth:
     run: True
     output_file_prefix: "mosdepth_results"
     bin_size: 1000
     mode: "fast"
     region_name: "LPA"
     work_dir: "path/to/mosdepth_work/"
     remove_intermediate: True

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Parameter
     - Description
   * - ``run``
     - Set ``True`` to run mosdepth.
   * - ``output_file_prefix``
     - Prefix for the aggregated depth output file.
   * - ``bin_size``
     - Genomic bin size in base pairs. ``1000`` (1 kb) is recommended —
       smaller bins increase resolution but memory and runtime.
   * - ``mode``
     - ``"fast"`` skips per-base depth (faster, less memory) or ``"full"``
       for per-base output. ``"fast"`` is sufficient for normalization.
   * - ``region_name``
     - Label applied to the VNTR region in output files (e.g. ``"LPA"``).
   * - ``work_dir``
     - Scratch directory for per-sample mosdepth intermediate files.
       Should be on fast local or scratch storage.
   * - ``remove_intermediate``
     - If ``True``, per-sample mosdepth files are deleted after aggregation
       to save disk space.

Step 3a — ``mosdepth.normalize``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Normalizes the depth matrix: filters low/high-coverage bins, applies a
repeat-mask, and z-scores across individuals. See
:doc:`algorithms/normalization` for the full method.

.. code-block:: yaml

   mosdepth:
     normalize:
       run: True
       min_depth: 20
       max_depth: 100
       top_frac: 0.1
       output_file_prefix: "mosdepth_results_normalized"
       repeat_mask_file: "path/to/repeat_mask.bed"

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Parameter
     - Description
   * - ``run``
     - Set ``True`` to normalize after mosdepth.
   * - ``min_depth``
     - Bins with median depth below this value are excluded.
   * - ``max_depth``
     - Bins with median depth above this value are excluded.
   * - ``top_frac``
     - Fraction of bins with the highest variance to retain for normalization
       (e.g. ``0.1`` keeps the top 10 % most variable bins).
   * - ``output_file_prefix``
     - Prefix for the normalized depth matrix output.
   * - ``repeat_mask_file``
     - BED file of repetitive regions to exclude before normalization.
       Download UCSC RepeatMasker tracks from the
       `UCSC Table Browser <https://genome.ucsc.edu/cgi-bin/hgTables>`_
       (group: Repeats, track: RepeatMasker).

Step 3b — ``mosdepth.neighbors``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Finds read depth-matched nearest neighbors in the normalized depth space
using scikit-learn ``NearestNeighbors``. See :doc:`algorithms/neighbors`.

.. code-block:: yaml

   mosdepth:
     neighbors:
       run: True
       output_file_prefix: "neighbor_coverage"
       num_neighbors: 5
       zmax: 2.0
       sigma2_max: 1000

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Parameter
     - Description
   * - ``run``
     - Set ``True`` to run neighbor finding.
   * - ``output_file_prefix``
     - Prefix for the neighbor assignment output file.
   * - ``num_neighbors``
     - Number of nearest neighbors to retain per sample. ``5`` is a
       conservative default; increase for small or heterogeneous cohorts.
   * - ``zmax``
     - Clip normalized depth values to :math:`\pm` ``zmax`` before distance computation.
       Reduces the influence of outlier bins. Typical range: ``1.5``–``3.0``.
   * - ``sigma2_max``
     - Maximum allowed neighbor variance. Samples whose nearest-neighbor
       set has variance above this threshold are flagged. Increase if many
       samples are being excluded in diverse cohorts.

Step 4 — ``compute_diploid_genotypes``
---------------------------------------

Estimates diploid copy number for each sample using the neighbor-normalized
depth ratio. For each sample, its locus depth is scaled by the mean depth of
its nearest neighbors at the same locus. See :doc:`algorithms/dipcn`.

.. code-block:: yaml

   compute_diploid_genotypes:
     run: True
     output_file_prefix: "diploid_genotypes"

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Parameter
     - Description
   * - ``run``
     - Set ``True`` to estimate diploid copy numbers.
   * - ``output_file_prefix``
     - Prefix for the diploid CN output file.

Step 5 — ``compute_haploid_genotypes``
---------------------------------------

Decomposes diploid CN estimates into haplotype-specific copy numbers using an
iterative IBD-based phasing algorithm. Requires a phased VCF and an IBS neighbors
file pre-computed by ``computeIBSpbwt``. See :doc:`algorithms/hi_inference` and
:doc:`ibs_ibd`.

.. code-block:: yaml

   compute_haploid_genotypes:
     run: True
     output_file_prefix: "haploid_genotypes"
     ibs_output: "path/to/ibs_neighbors_chr6.tsv.gz"
     min_neighbors: 1
     max_neighbors: 10
     n_iters: 100

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Parameter
     - Description
   * - ``run``
     - Set ``True`` to run haplotype inference.
   * - ``output_file_prefix``
     - Prefix for the haploid CN output files.
   * - ``ibs_output``
     - Path to the IBS/IBD neighbors file produced by ``computeIBSpbwt``.
       See :doc:`ibs_ibd` for how to generate this file.
   * - ``min_neighbors``
     - Minimum number of IBS neighbors required to include a sample.
       Samples with fewer neighbors are excluded from phasing.
   * - ``max_neighbors``
     - Maximum number of IBS neighbors used per haplotype per iteration.
       Controls the breadth of the phasing signal.
   * - ``n_iters``
     - Number of iterations for the haplotype update loop. Convergence is
       typically reached within 20–50 iterations; ``100`` is a safe default.

.. note::

   ``ibs_output`` is the only required external file for this step. It must be
   generated before running the pipeline. See :doc:`ibs_ibd` for full instructions
   and :doc:`examples/IBS_example` for an annotated example script.

----

Full Example Config
-------------------

.. literalinclude:: ../../grid/example_config.yaml
   :language: yaml
   :linenos:
   :caption: grid/example_config.yaml
