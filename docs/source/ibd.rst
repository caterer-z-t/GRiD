Computing IBD Neighbors with iLASH
===================================

Step 7 of GRiD (haplotype inference) can use identity-by-descent (IBD) haplotype
neighbors computed by `iLASH <https://github.com/roohy/iLASH>`_, a scalable
local alignment-based IBD detection tool. IBD segments represent true
identity-by-descent rather than identity-by-state, making them more specific
haplotype neighbors for phasing — particularly for individuals who lack long
IBS matches.

For IBS-based neighbors using ``computeIBSpbwt``, see :doc:`ibs`.

.. note::

   iLASH is not included in the GRiD repository and must be available on your
   ``PATH`` (e.g. via ``module load ilash`` on supported HPC clusters).

----

Overview
--------

The iLASH workflow for GRiD requires three stages before running iLASH itself:

1. **Phase your VCF** — phased genotypes are required (e.g. from
   `GnoMix <https://github.com/AI-sandbox/gnomix>`_ or a reference panel).
2. **Convert VCF → PED/MAP** — iLASH takes PLINK text-format (``.ped/.map``)
   input, not VCF or BGEN. Use ``vcf2ped_zc.py`` to preserve phase information.
3. **Add genetic map coordinates** — iLASH needs cM positions. Use
   ``add_gen_map.py`` to annotate the MAP file with a genetic map.

----

Input Requirements
------------------

**Phased VCF** — iLASH requires phased genotypes. If your VCF is unphased,
phase it first with a tool such as GnoMix, Eagle, or SHAPEIT. The phased VCF
must be indexed and cover the chromosome of interest.

**Genetic map** — a tab-separated file with physical and cM positions for
the target chromosome. Download pre-built hg38 maps:

.. code-block:: bash

   wget https://alkesgroup.broadinstitute.org/Eagle/downloads/tables/genetic_map_hg38.tar.gz
   tar -xzf genetic_map_hg38.tar.gz

----

Stage 1 — Convert VCF to PED/MAP
----------------------------------

iLASH reads PLINK text format (``.ped/.map``). The ``vcf2ped_zc.py`` helper
converts a phased VCF while preserving haplotype phase. It lives in
``grid/utils/helper_dir/``:

.. code-block:: bash

   # Count variants in the VCF
   SNP_COUNT=$(bcftools view -H phased_chr6.vcf.gz | wc -l)

   python grid/utils/helper_dir/vcf2ped_zc.py \
       --vcf       phased_chr6.vcf.gz \
       --out       gnomix_phased_chr6 \
       --snp-count $SNP_COUNT

This produces ``gnomix_phased_chr6.ped`` and ``gnomix_phased_chr6.map``.

.. note::

   Use ``vcf2ped_zc.py`` rather than standard PLINK conversion tools —
   standard converters do not guarantee phase order is preserved in the
   PED file, which is required for iLASH to assign haplotypes correctly.

----

Stage 2 — Add Genetic Map Coordinates
----------------------------------------

iLASH uses cM positions from the MAP file to compute segment lengths.
Annotate the MAP file produced in Stage 1 using ``add_gen_map.py`` from
``grid/utils/helper_dir/``:

.. code-block:: bash

   python grid/utils/helper_dir/add_gen_map.py \
       --map         gnomix_phased_chr6.map \
       --genetic-map genetic_map_hg38_chr6.txt.gz \
       --out         gnomix_phased_chr6_mapped

This produces ``gnomix_phased_chr6_mapped.map`` with cM coordinates filled in.

----

Stage 3 — iLASH Configuration File
-------------------------------------

iLASH is configured via a plain-text parameter file. Key parameters for GRiD use:

.. code-block:: text

   # iLASH parameter file for LPA KIV-2 (chr6)
   map  /path/to/gnomix_phased_chr6_mapped.map
   ped  /path/to/gnomix_phased_chr6.ped
   output /path/to/output.match

   # SNP-based slicing
   auto_slice 1
   slice_size 350
   step_size 350
   min_length 0.5
   slice_length 0.5
   cm_overlap 1

   # LSH parameters
   shingle_size 15
   shingle_overlap 0
   perm_count 50
   bucket_count 10
   interest_threshold 0.50
   match_threshold 0.70
   minhash_threshold 55

   # Threading
   max_thread 3

.. list-table::
   :header-rows: 1
   :widths: 28 72

   * - Parameter
     - Description
   * - ``map`` / ``ped``
     - Paths to the genetic-map-annotated MAP file and phased PED file.
   * - ``output``
     - Path for the output ``.match`` file.
   * - ``auto_slice``
     - Set ``1`` to use SNP-count-based window slicing (recommended).
   * - ``slice_size`` / ``step_size``
     - Window size and step in number of SNPs.
   * - ``min_length`` / ``slice_length``
     - Minimum IBD segment length (cM) and slice length (cM).
   * - ``shingle_size`` / ``perm_count`` / ``bucket_count``
     - LSH parameters controlling sensitivity vs. speed.
       Increasing ``perm_count`` improves recall at the cost of runtime.
   * - ``interest_threshold`` / ``match_threshold``
     - Candidate and final match score thresholds written to the output file.
       Set permissively here (e.g. ``0.70``) and let GRiD filter at load time.
   * - ``minhash_threshold``
     - Minimum number of MinHash collisions to consider a pair. Lower values
       increase sensitivity but also false-positive rate.
   * - ``max_thread``
     - Number of CPU threads. Match to your SLURM ``--cpus-per-task``.

.. tip::

   Run iLASH with a permissive ``match_threshold`` (e.g. ``0.70``) and let
   GRiD filter at load time via ``min_match`` and ``match_at_min_length``.
   This avoids re-running iLASH when experimenting with different quality
   thresholds during phasing.

----

Running iLASH
--------------

Once the config file is written, run iLASH:

.. code-block:: bash

   ilash /path/to/config

Full SLURM submission script:

.. code-block:: bash

   #!/bin/bash
   #SBATCH --job-name=ilash
   #SBATCH --output=slurm/%x.out
   #SBATCH --error=slurm/%x.err
   #SBATCH --time=4:00:00
   #SBATCH --mem=300G
   #SBATCH --cpus-per-task=3

   set -euo pipefail

   module load bcftools ilash anaconda
   conda activate gnomix

   OUTDIR="path/to/output"
   STUDY="my_cohort"
   SCRIPTS="path/to/grid/scripts"

   mkdir -p "$OUTDIR/$STUDY/data"

   # Decompress phased VCF if needed
   PHASED_VCF="$OUTDIR/$STUDY/data/gnomix_output/query_file_phased.vcf.gz"
   gunzip -c "$PHASED_VCF" > "${PHASED_VCF%.gz}"
   PHASED_VCF="${PHASED_VCF%.gz}"

   SNP_COUNT=$(bcftools view -H "$PHASED_VCF" | wc -l)

   # Stage 1: VCF → PED/MAP
   python "$SCRIPTS/grid/utils/helper_dir/vcf2ped_zc.py" \
       --vcf "$PHASED_VCF" \
       --out "$OUTDIR/$STUDY/data/gnomix_phased_${STUDY}" \
       --snp-count "$SNP_COUNT"

   # Stage 2: add genetic map coordinates
   python "$SCRIPTS/grid/utils/helper_dir/add_gen_map.py" \
       --map         "$OUTDIR/$STUDY/data/gnomix_phased_${STUDY}.map" \
       --genetic-map "$SCRIPTS/files/genetic_map_hg38_chr6.txt.gz" \
       --out         "$OUTDIR/$STUDY/data/gnomix_phased_${STUDY}_mapped"

   # Stage 3: write iLASH config and run
   cat > "$OUTDIR/$STUDY/config" << EOT
   map    $OUTDIR/$STUDY/data/gnomix_phased_${STUDY}_mapped.map
   ped    $OUTDIR/$STUDY/data/gnomix_phased_${STUDY}.ped
   output $OUTDIR/$STUDY/data/output.match
   auto_slice 1
   slice_size 350
   step_size 350
   min_length 0.5
   slice_length 0.5
   cm_overlap 1
   shingle_size 15
   shingle_overlap 0
   perm_count 50
   bucket_count 10
   interest_threshold 0.50
   match_threshold 0.70
   minhash_threshold 55
   max_thread $SLURM_CPUS_PER_TASK
   EOT

   ilash "$OUTDIR/$STUDY/config"

   # Clean up intermediate files
   rm "$OUTDIR/$STUDY/data/gnomix_phased_${STUDY}.map" \
      "$OUTDIR/$STUDY/data/gnomix_phased_${STUDY}.ped" \
      "$OUTDIR/$STUDY/data/gnomix_phased_${STUDY}_mapped.map"

.. note::

   Memory usage scales with cohort size and chromosome length.
   For a 2,500-sample cohort on chr6, expect 200–300 GB RAM and 2–4 hours
   at 3 threads.

----

Output Format
-------------

iLASH writes a tab-separated ``.match`` file with one IBD segment per row
and no header line:

.. code-block:: text

   NA19919  NA19919_1  NA21144  NA21144_1  chr6  166290821  166437824  6:166290821:C:A  6:166437824:T:C  0.541199  0.926214
   HG03563  HG03563_0  NA21144  NA21144_0  chr6  160157539  160706689  6:160157539:G:T  6:160706689:T:C  0.515263  0.798077

.. list-table::
   :header-rows: 1
   :widths: 8 20 72

   * - Col
     - Name
     - Description
   * - 1
     - ``ID1``
     - Sample identifier for individual 1
   * - 2
     - ``HAP_ID1``
     - Haplotype identifier: ``{ID1}_{hap}`` where hap is ``0`` or ``1`` (0-indexed)
   * - 3
     - ``ID2``
     - Sample identifier for individual 2
   * - 4
     - ``HAP_ID2``
     - Haplotype identifier: ``{ID2}_{hap}`` where hap is ``0`` or ``1`` (0-indexed)
   * - 5
     - ``chrom``
     - Chromosome (e.g. ``chr6``)
   * - 6
     - ``bp_start``
     - Start position of the IBD segment (bp, same build as input)
   * - 7
     - ``bp_end``
     - End position of the IBD segment (bp)
   * - 8
     - ``snp_start``
     - Variant ID at the start of the segment
   * - 9
     - ``snp_end``
     - Variant ID at the end of the segment
   * - 10
     - ``length_cM``
     - Segment length in centiMorgans
   * - 11
     - ``match_score``
     - iLASH match score in [0, 1]. The default IBD calling threshold in iLASH
       is ≥ 0.99, but lower values are retained when ``match_threshold`` /
       ``interest_threshold`` are relaxed.

.. note::

   iLASH haplotype indices are **0-indexed** (``_0`` and ``_1``), whereas
   ``computeIBSpbwt`` uses **1-indexed** haplotypes (``hap=1`` or ``hap=2``).
   GRiD's ``hi_inference`` handles both formats automatically.

----

Passing the Output to GRiD
---------------------------

Point to the ``.match`` file in your config under ``ibd_output``:

.. code-block:: yaml

   compute_haploid_genotypes:
     run: True
     output_file_prefix: "haploid_genotypes"
     method: "ibd"
     phasing_method: "graph_prop"
     ibd_output: "path/to/output.match"
     alpha: 0.7
     n_iters: 20
     max_neighbors: 20
     min_length: 2.0
     min_match: 0.95

The top-level ``start_bp`` and ``end_bp`` config keys (the VNTR coordinates)
are used to assign IBD segments to Tier-1 (overlapping the VNTR) or Tier-2
(elsewhere on the chromosome) priority buckets when loading neighbors.
See :ref:`tiered-neighbors` for details on the two-tier system and
:ref:`adaptive-threshold` for the length-adaptive match threshold.

For a full description of all IBD-related config parameters and phasing
algorithms, see :doc:`algorithms/hi_inference`.