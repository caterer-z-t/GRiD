---
title: 'Genomic Repeat inference from Depth (GRiD): A Pipeline for Haplotype-Resolved Large VNTR Copy Number Estimation'
tags:
  - Python
  - Variable Nucleotide Tandem Repeat (VNTR)
  - Copy Number Variation (CNV)
  - Identity By Descent (IBD)
  - Haplotype Inference
  - Lipoprotein(a)
  - Cardiovascular Disease
authors:
  - name: Zachary Caterer
    orcid: 0000-0001-9019-0730
    corresponding: true
    affiliation: "1, 2"
  - name: Meng Lin
    orcid: 0000-0003-4603-0718
    affiliation: 1
  - name: Qiang An
    affiliation: 3
  - name: Maisy Brasher
    orcid: 0000-0002-4020-6551
    affiliation: 1
  - name: Joanne Cole
    orcid: 0000-0001-9520-2788
    affiliation: 1
  - name: Ethan Lange
    affiliation: 1
    orcid: 0000-0001-7075-4287
  - name: Christopher Gignoux
    affiliation: 1
  - name: Mariaelisa Graff
    affiliation: 3
    orcid: 0000-0001-6380-1735
  - name: Christy Avery
    affiliation: 3
    orcid: 0000-0002-1044-8162
  - name: Maggie Stanislawski
    affiliation: 1
    orcid: 0000-0001-7768-8868
    corresponding: true
affiliations:
 - name: Department of Biomedical Informatics, University of Colorado Anschutz, Aurora, CO USA
   index: 1
   ror: "03wmf1y16"
 - name: Department of Chemical and Biological Engineering, University of Colorado Boulder, Boulder, CO USA
   index: 2
   ror: "02ttsq026"
 - name: Department of Epidemiology, University of North Carolina Chapel Hill, Chapel Hill, NC USA
   index: 3
   ror: "0130frc33"
date: \today
bibliography: paper.bib
---

# Summary

Variable number tandem repeats (VNTRs) are genomic regions in which a sequence motif is repeated in tandem a variable number of times across individuals. Large VNTRs — those with repeat units spanning hundreds to thousands of base pairs — are a substantial source of structural variation in the human genome and have been linked to complex traits and disease [@mukamel2021]. Because their repeat units far exceed the length of standard short sequencing reads, large VNTRs cannot be genotyped by conventional variant-calling approaches, and haplotype-resolved copy number estimation from short-read whole-genome sequencing (WGS) remains a methodological gap. @mukamel2021 identified 734 protein-coding VNTRs with elevated identity-by-descent (IBD) rates across the human genome, many of which are likely functional but lack accessible copy number estimation tools for cohort-scale short-read data.

A particularly well-characterized example is the Kringle IV type-2 (KIV-2) VNTR within the *LPA* gene, which is the primary genetic determinant of Lipoprotein(a) [Lp(a)] concentration [@utermann1987; @mukamel2021]. Lp(a) is a causal, independent, and highly heritable ($\text{h}^2 >70\%$) cardiovascular disease (CVD) risk factor elevated in an estimated 1.5 billion people globally, with concentrations differing substantially across ancestry groups [@tsimikas2017; @kronenberg2022; @nordestgaard2010]. The KIV-2 repeat unit spans ~5.5 kb, individuals carry 1–40 copies per haplotype, and copy number is inversely correlated with Lp(a) concentration, making accurate haplotype-resolved KIV-2 CNV estimation central to ancestry-aware polygenic risk scores and causal inference studies [@utermann1987; @clarke2009; @kamstrup2010; @mukamel2021].

**GRiD** (**G**enomic **R**epeat **i**nference from **D**epth) is an open-source Python pipeline for haplotype-resolved large VNTR copy number estimation from short-read WGS data. GRiD combines read-depth normalization, read-depth profile-matched nearest-neighbor estimation, and IBD-based haplotype inference [@hujoel2026] into a single installable, multi-threaded, YAML-configured workflow. GRiD was developed and validated for the *LPA* KIV-2 locus, but its core framework is locus-agnostic: any large VNTR in a cohort with sufficient sample size can be targeted by supplying the relevant genomic coordinates in the configuration file. The pipeline is designed to be accessible to researchers without extensive bioinformatics expertise.

# Statement of Need

Large VNTRs present a shared set of computational challenges that place them outside the scope of existing repeat genotyping software. Their repeat units are orders of magnitude larger than the short tandem repeats (STRs, 1–6 bp motifs) targeted by most tools, and in high-copy-number individuals the VNTR region can span hundreds of kilobases, making read-pair spanning approaches infeasible with standard 150 bp paired-end reads. Haplotype resolution adds a further challenge: phasing cannot rely on long reads — which remain unavailable for most large epidemiological cohorts — or on pre-existing phased variant data at the locus. GRiD is designed to address these challenges directly through three features: (1) a read-depth normalization and nearest-neighbor framework that is locus-agnostic and parameterized entirely through a configuration file, (2) an IBD-based iterative phasing algorithm [@hujoel2026] that decomposes diploid copy number into haplotype-specific estimates using cohort relatedness alone, and (3) a single installable pipeline that encapsulates the full workflow from short-read WGS to haplotype-resolved copy number in a reproducible, multi-threaded form.

For the specific case of *LPA* KIV-2 — GRiD's primary development and validation target — several tools have addressed parts of this problem, but none provides the complete short-read-to-haplotype-resolved workflow in open-source, cohort-scalable form. KILDA [@molitor2025] estimates diploid KIV-2 copy number from FASTQ files using an alignment-free k-mer approach but does not resolve haplotype-specific copy numbers. The Illumina DRAGEN LPA Caller provides KIV-2 estimation but requires the proprietary DRAGEN platform. vntr-calling-nf [@vntrcallingnf] offers a Nextflow-based pipeline for VNTR variant calling but targets coding variant detection rather than copy number quantification and has been primarily validated on whole-exome sequencing data. As a result, researchers studying Lp(a) at cohort scale with short-read WGS have largely relied on in-house scripts or unpublished internal pipelines, creating barriers to reproducibility and cross-cohort comparability.

# State of the Field

Several mature tools exist for genotyping tandem repeats from short-read sequencing data. ExpansionHunter [@dolzhenko2019] and GangSTR [@mousavi2021] use probabilistic models of read-pair orientation and coverage to estimate repeat lengths, but are designed for STRs and short VNTRs where reads can partially or fully span the repeat unit. HipSTR [@willems2017] performs haplotype-level STR genotyping using phased reads but is similarly constrained to short repeat units accessible by individual reads. Tools such as TRTools [@mousavi2021trtools] provide downstream analysis utilities for STR calls but do not extend to large-VNTR loci. adVNTR [@bakhtiari2018] uses hidden Markov models to genotype VNTRs genome-wide and supports both Illumina and PacBio data, but is designed for VNTRs of up to a few hundred base pairs where reads can overlap the repeat boundary, which is qualitatively different from the 5.6 kb KIV-2 repeat unit. More recently, long-read platforms have enabled direct assembly of complex repeat regions [@didericksen2024], and Kivvi [@pacbio2026kivvi] specifically targets KIV-2 from PacBio HiFi data; however, long-read data are not yet available for most large epidemiological cohorts, which remain dominated by short-read WGS.

For the *LPA* KIV-2 locus, KILDA [@molitor2025] provides an alignment-free diploid copy number estimate from short-read FASTQ files using a k-mer counting strategy calibrated against long-read assemblies, demonstrating strong concordance with assembly-based diploid estimates. However, no existing tool addresses haplotype-level phasing for KIV-2 or for large VNTRs more broadly from short-read data. The foundational algorithmic framework underlying GRiD was established by @mukamel2021, who demonstrated that population-normalized read depth is the most tractable approach for large VNTR copy number estimation from short-read cohort data, and by @hujoel2026, who developed an IBD-based iterative algorithm for decomposing diploid copy number into haplotype-specific contributions without requiring long reads or pre-existing phased variant data. GRiD integrates and operationalizes both methodological contributions into a single accessible software package with read-depth profile-aware nearest-neighbor normalization, a full configuration system, and multi-cohort scalability.

# Software Design

GRiD is implemented in Python (≥3.8) and follows a modular architecture organized into seven sequential pipeline steps, each implemented as an independent utility module under `grid/utils/`. The pipeline is orchestrated by `grid/pipeline.py` and driven by a user-supplied YAML configuration file, which specifies file paths, genomic coordinates, per-step parameters, and which steps to execute. The command-line interface is implemented with Click [@click] and provides an entry point: `grid wgs` for whole-genome sequencing data.

The Whole Genome Sequencing pipeline proceeds as follows:

1. **Index verification / creation**: Ensures CRAM/BAM index files (.crai/.bai) exist for all input samples, creating them with `samtools` [@li2009] if needed.
2. **Read counting**: Counts properly paired reads in the target VNTR region using pysam [@pysam], filtering by mapping quality and SAM flag.
3. **Coverage estimation**: Runs mosdepth [@pedersen2018] across the genome in binned mode to produce per-sample depth profiles.
4. **Coverage normalization**: Normalizes coverage within individuals (by sample mean depth) and across individuals (by a z-score transformation), then filters to high-variance regions enriched for VNTR signal.
5. **Nearest-neighbor identification**: Uses scikit-learn [@pedregosa2011] to compute Euclidean distances in normalized depth space and identify the top-N genomically similar neighbors per individual, providing a cohort-matched reference for CNV normalization.
6. **Diploid copy number estimation**: Computes per-exon diploid copy number for each individual by normalizing their read counts relative to their nearest neighbors' read counts, accounting for sample-specific sequencing depth.
7. **Haplotype inference**: Applies an IBD-based iterative phasing algorithm [@hujoel2026] that uses haplotype-matched neighbors to decompose diploid copy number estimates into haplotype-specific contributions.

All multi-sample steps are parallelized using Python's `concurrent.futures.ThreadPoolExecutor`, and progress is reported via the Rich terminal library [@rich]. Configuration validation at pipeline startup provides early, informative error messages before any compute-intensive steps run.

# Research Impact Statement

GRiD has been validated across multiple large-scale, multi-ancestry whole-genome sequencing cohorts representing tens of thousands of individuals spanning Hispanic/Latino, African American, and European American ancestry groups. Its nearest-neighbor normalization approach naturally accounts for ancestry and sequencing batch effects without requiring explicit ancestry labels, making it well-suited for population-diverse studies where CNV tools calibrated on European-ancestry reference panels may systematically underperform.

GRiD has attracted national and international interest from research groups across multiple conferences, with inquiries from groups seeking to apply large VNTR copy number estimation to their own short-read WGS cohorts. It is designed for accessibility, with comprehensive documentation of configuration, outputs, and algorithmic details to lower the barrier to entry for epidemiologists and clinicians without extensive bioinformatics training.

# Mathematics

**Coverage estimation.** Weighted mean depth across the VNTR region from mosdepth bins is computed as:

$$
\text{Coverage} =
\left\lfloor
100 \cdot
\frac{
\sum_{i \in R} \left( \bar{C}_i \cdot \bigl( \min(e, r_{i,e}) - \max(s, r_{i,s}) \bigr) \right)
}{
\sum_{i \in R} \bigl( \min(e, r_{i,e}) - \max(s, r_{i,s}) \bigr)
}
\right\rceil
$$

where $R$ is the set of bins overlapping the region $[s, e]$, $\bar{C}_i$ is the mean depth of bin $i$, and $r_{i,s}$, $r_{i,e}$ are the bin boundaries.

**Normalization.** Let $D \in \mathbb{R}^{N \times M}$ be the depth matrix for $N$ individuals across $M$ genomic bins. Within-individual normalization gives $D^{(1)}_{ij} = D_{ij} / \bar{D}_i$, where $\bar{D}_i$ is the individual mean. Across-individual normalization then yields:

$$
D^{\mathrm{norm}}_{ij} = \frac{D^{(1)}_{ij} - \mu_j}{\sqrt{\mu_j}}
$$

where $\mu_j$ is the population mean of $D^{(1)}_{\cdot j}$. Regions are filtered by a variance ratio $\sigma_j^2 / \mu_j$, and the top fraction of high-variance regions are retained for neighbor computation.

**Neighbor-normalized diploid copy number.** For individual $i$ with read count $r_i$ and depth scale $s_i$, the diploid copy number estimate is:

$$
\widehat{\text{dipCN}}_i = \frac{r_i / s_i}{\frac{1}{|N_i|} \sum_{j \in N_i} r_j / s_j}
$$

where $N_i$ is the set of nearest neighbors identified from the normalized depth space.

**Haplotype inference.** Following @hujoel2026, haplotype-specific copy numbers $h_{i,1}$ and $h_{i,2}$ are iteratively estimated by updating each haplotype's value proportionally to the mean copy number of IBD-matched haplotype neighbors:

$$
h_{i,k}^{(t+1)} = \widehat{\text{dipCN}}_i \cdot \frac{\bar{h}_{N_{i,k}}^{(t)}}{\bar{h}_{N_{i,1}}^{(t)} + \bar{h}_{N_{i,2}}^{(t)}}
$$

where $\bar{h}_{N_{i,k}}^{(t)}$ is the mean haplotype copy number of IBD neighbors on haplotype $k$ at iteration $t$.

# AI Usage Disclosure

GitHub Copilot and Anthropic Claude were used to assist with code development and documentation writing during the preparation of this software and manuscript. All AI-generated content was reviewed line-by-line by the authors, with corrections and edits applied as necessary to ensure accuracy, correctness, and consistency with the underlying methodology. No AI-generated code or text was incorporated without direct author verification.

# Acknowledgements

This work was supported by [GRANT INFORMATION].

# References
