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
    affiliation: "1, 2, 3"
  - name: Meng Lin 
    orcid: 0000-0003-4603-0718 
    affiliation: 1
  - name: Qiang An 
    affiliation: 4
  - name: Maizy S. Brasher
    orcid: 0000-0002-4020-6551
    affiliation: 1
  - name: Elbay Aliyev
    orcid: 0000-0002-6469-1854
    affiliation: 1
  - name: Harriet Dashnow
    orcid: 0000-0001-8433-6270
    affiliation: 1
  - name: Joanne B. Cole
    orcid: 0000-0001-9520-2788
    affiliation: "1, 3"
  - name: Ethan Lange 
    affiliation: 1 
    orcid: 0000-0001-7075-4287
  - name: Mariaelisa Graff 
    affiliation: 4 
    orcid: 0000-0001-6380-1735
  - name: Christy L. Avery 
    affiliation: 4 
    orcid: 0000-0002-1044-8162
  - name: Christopher R. Gignoux 
    affiliation: "1, 3"
    orcid: 0000-0001-9728-6567
  - name: Maggie Stanislawski 
    affiliation: "1, 3"
    corresponding: true

affiliations:
 - name: Department of Biomedical Informatics, University of Colorado Anschutz, Aurora, CO USA
   index: 1
   ror: "03wmf1y16"
 - name: Department of Chemical and Biological Engineering, University of Colorado Boulder, Boulder, CO USA
   index: 2
   ror: "02ttsq026"
 - name: Interdisciplinary Quantitative Biology PhD Program, Biofrontiers Institute, University of Colorado Boulder, Boulder, CO USA 
   index: 3 
   ror: "02ttsq026"
 - name: Department of Epidemiology, University of North Carolina Chapel Hill, Chapel Hill, NC USA
   index: 4
   ror: "0130frc33"
date: \today
bibliography: paper.bib
---

# Summary

Variable number tandem repeats (VNTRs) are genomic regions in which a sequence motif is repeated in tandem a variable number of times across individuals. Large VNTRs, repeat units spanning hundreds to thousands of base pairs, are a substantial source of structural variation in the human genome and have been linked to complex traits and disease [@mukamel2021]. Because repeat units in these VNTRs can be more than 35 times the length of the standard 150-bp short sequencing reads, large VNTRs cannot be genotyped by conventional variant-calling approaches. Furthermore, accurate copy number estimation from short-read whole-genome sequencing (WGS) remains challenging. A well-characterized example is the Kringle IV type-2 (KIV-2) VNTR within the LPA gene, which is the primary genetic determinant of lipoprotein(a) [Lp(a)] concentration [@utermann1987; @mukamel2021]. Lp(a) is a causal and highly heritable ($\text{h}^2 >70\%$) atherosclerotic cardiovascular disease (ASCVD) risk factor elevated in an estimated 1.5 billion people globally [@reyes2022lipoprotein]. Both Lp(a) concentrations and effect of LPA genetic variants vary substantially across ancestry groups, underscoring the importance of ancestry-aware genetic analyses and risk prediction [@tsimikas2017; @kronenberg2022; @nordestgaard2010]. The KIV-2 repeat unit spans ~5.5 kb, with individuals carrying 1–40 copies per haplotype, and KIV-2 copy number is inversely correlated with Lp(a) concentration [@utermann1987; @clarke2009; @kamstrup2010; @mukamel2021]. Accurate estimation of KIV-2 copy number therefore has important implications for genetic studies of Lp(a) and cardiovascular disease, including analyses of ancestry-related variation and genetic risk [@mukamel2021; @kronenberg2022]. While this VNTR serves as a region of great interest, prior work has been limited to custom pipelines, limiting applications by the broader field.

**GRiD** (**G**enomic **R**epeat **i**nference from **D**epth) is an open-source Python pipeline that addresses this methodological gap. GRiD estimates diploid and haplotype-resolved VNTR copy number from short-read WGS data. GRiD combines read-depth normalization, read depth profile-matched nearest-neighbor estimation, and identity-by-descent (IBD)-based haplotype inference [@hujoel2026]. Although GRiD was developed and validated for the LPA KIV-2 locus, its core framework is locus-agnostic and applicable to any large VNTR with sufficient sample size. The pipeline is designed to be accessible to researchers without extensive bioinformatics expertise.

# Statement of need

Large VNTRs present a shared set of computational challenges that place them outside the scope of existing repeat genotyping software and standard genomic analysis pipelines. Their repeat units are orders of magnitude larger than the short tandem repeats (STRs) targeted by most available tools [@dolzhenko2019; @mousavi2021; @willems2017; @bakhtiari2018], and in individuals with high-copy-number alleles the VNTR region can span hundreds of kilobases, making read-pair spanning approaches difficult with standard paired-end reads. Haplotype resolution from short-read data can be important for downstream genomic analyses that depend allele-specific copy number, including analyses of statistical associations, genomic coverage, and haplotype based tract, particularly given the expense and difficulty of obtaining long-read data for large epidemiological cohorts. GRiD is designed to address these challenges directly through three features: (1) a read-depth normalization and nearest-neighbor framework that is locus-agnostic and parameterized entirely through a configuration file, (2) an IBD-based iterative phasing algorithm [@hujoel2026] that decomposes diploid copy number into haplotype-specific estimates using cohort relatedness, and (3) a single installable pipeline that encapsulates the full workflow from short-read WGS to haplotype-resolved copy number in an intuitive, reproducible, multi-threaded form. 

Although several software packages address aspects of tandem repeat analysis, no existing open-source tool provides an end-to-end workflow for estimating both diploid and haplotype-resolved copy number from short-read whole-genome sequencing data at modern biobank scales. Consequently, many studies have relied on custom analysis pipelines, limiting reproducibility, portability, and cross-cohort compatibility. GRiD fills this gap by providing a configurable, documented, and openly available implementation of these methods. 

# State of the field

Several mature software packages exist for genotyping tandem repeats from short-read sequencing data. ExpansionHunter [@dolzhenko2019], and HipSTR [@willems2017] were developed primarily for short tandem repeats (STRs) whereas GangSTR [@mousavi2021] was developed for short VNTRs (defined by a repeat motif length of up to 20bp). However, these methods are designed for repeat lengths where short reads can still provide sufficient direct evidence about repeat structure. Likewise, TRTools [@mousavi2021trtools] provides downstream analysis utilities for STR genotypes but does not extend to large VNTR copy number estimation. adVNTR [@bakhtiari2018] supports genome-wide VNTR genotyping using hidden Markov models but is designed for repeat units that are substantially shorter (<150 bp) than loci such as the 5.6 kb LPA KIV-2 repeat. Earlier methods such as VNTRseek [@VNTRseek] detected shorter VNTRs by comparing read-derived repeat counts to a reference catalog, while GtTR [@GtTR] estimated absolute target VNTR copy number by scaling short-read depth across locus boundaries against a baseline genotype derived from long reads, though it remains unable to phase individual diploid alleles. More recently, danbing-tk [@lu2021profiling] introduced repeat-pangenome graphs for genome-wide profiling of VNTR length and motif composition from short-read sequencing data. This approach represents VNTR sequence diversity using haplotype-resolved assemblies and estimates VNTR dosage from graph-specific k-mer counts. Recent advances in long-read sequencing have enabled direct characterization of complex repeat regions [@didericksen2024]. For example, Kivvi [@pacbio2026kivvi] estimates KIV-2 copy number from PacBio HiFi sequencing data. However, most large population cohorts (eg, 1000G [@auton2015global], TOPMed [@taliun2021sequencing], and the Genotype-Tissue Expression (GTEx) project [@GTConsortium2017nature]) remain based on short-read whole-genome sequencing, motivating methods that can estimate VNTR copy number without requiring long-read data. KILDA [@molitor2025] estimates KIV-2 diploid copy number directly from short-read FASTQ files using an alignment-free k-mer-counting approach normalized against non-repetitive sequences. The proprietary Illumina DRAGEN LPA Caller similarly estimates KIV-2 copy number but requires the FPGA-accelerated DRAGEN platform. vntr-calling-nf [@vntrcallingnf] provides a Nextflow workflow for VNTR analysis but focuses on coding variant detection rather than large-VNTR copy number quantification and has primarily been validated using whole-exome sequencing data. Collectively, these tools do not provide an end-to-end solution for estimating haplotype-resolved copy number of large VNTRs from short-read whole-genome sequencing data at population scale.

# Software design

GRiD integrates approaches from prior work for estimating large VNTR copy number from short-read sequencing [@mukamel2021; @hujoel2026]. It provides a configurable, end-to-end workflow for estimating both diploid and haplotype-resolved large VNTR copy number from short-read whole-genome sequencing data at cohort scale. The pipeline is designed to be accessible to researchers without extensive bioinformatics expertise. 

GRiD is implemented in Python ($\geq$ 3.8) and follows a modular architecture organized into seven sequential pipeline steps, each implemented as an independent utility module under `grid/utils/`. The pipeline is orchestrated by `grid/pipeline.py` and driven by a user-supplied YAML configuration file, which specifies file paths, genomic coordinates, per-step parameters, and which steps to execute. The command-line interface is implemented with Click [@click] and provides an entry point: `grid wgs` for whole-genome sequencing data.

The Whole Genome Sequencing pipeline proceeds as follows:

1. **Index verification / creation**: Ensures CRAM/BAM index files (`.crai`/`.bai`) exist for all input samples, creating them with `samtools` [@li2009] if needed.
2. **Read counting**: Counts properly paired reads in the target VNTR region using pysam [@pysam], filtering by mapping quality and SAM flags.
3. **Coverage estimation**: Runs mosdepth [@pedersen2018] across the genome in binned mode to produce per-sample depth profiles.
4. **Coverage normalization**: Normalizes coverage within individuals (by sample mean depth) and across individuals (by a z-score transformation), then filters to high-variance regions enriched for VNTR signal.
5. **Nearest-neighbor identification**: Uses scikit-learn [@pedregosa2011] to compute Euclidean distances in normalized depth space and identify the top-N genomically similar neighbors per individual, providing a cohort-matched reference for CNV normalization.
6. **Diploid copy number estimation**: Computes per-exon diploid copy number for each individual by normalizing their read counts relative to their nearest neighbors' read counts, accounting for sample-specific sequencing depth.
7. **Haplotype inference**: Applies an IBD-based iterative phasing algorithm [@hujoel2026] that uses haplotype-matched neighbors to decompose diploid copy number estimates into haplotype-specific contributions.

All multi-sample steps are parallelized using Python's `concurrent.futures.ThreadPoolExecutor`, and progress is reported via the Rich terminal library [@rich]. Configuration validation at pipeline startup provides early, informative error messages before any compute-intensive steps run.

# Research impact statement

GRiD has been applied to >50,000 study participants with whole-genome sequencing data from multiple TOPMed cohorts and 1000G cohorts, comprising tens of thousands of participants across African, European, Hispanic/Latino, and admixed American ancestry groups. These analyses form the basis of ongoing research manuscripts and demonstrate the applicability of GRiD to population-scale analyses and heterogeneous sequencing datasets.

GRiD provides a reproducible, open-source implementation of a previously custom analysis workflow, enabling large-VNTR copy number estimation from existing short-read whole-genome sequencing datasets. By packaging the complete workflow with documentation, example datasets, and automated installation, GRiD reduces the technical barriers to incorporating large VNTR analyses into population-scale genomic studies.

# Mathematics

The mathematical framework implemented in GRiD follows previously published methods. Coverage estimation, normalization, and diploid copy-number estimation follow Mukamel et al. [@mukamel2021], while haplotype-specific copy-number inference follows Hujoel et al. [@hujoel2026]. The equations below describe their implementation within GRiD; detailed methodological derivations and justification are provided in the original publications.

**Coverage estimation.** Following Mukamel et al. [@mukamel2021], weighted mean depth across the VNTR region from mosdepth [@mosdepth] bins is computed as:

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

where $R$ is the set of bins overlapping the region $[s,e]$, $\bar{C}_i$ is the mean depth of bin $i$, and $r_{i,s}$ and $r_{i,e}$ are the bin boundaries.

**Normalization.** Following Mukamel et al. [@mukamel2021], let $D \in \mathbb{R}^{N \times M}$ be the depth matrix for $N$ individuals across $M$ genomic bins. Within-individual normalization gives $D^{(1)}_{ij} = D_{ij} / \bar{D}_i$, where $\bar{D}_i$ is the individual mean. Across-individual normalization then yields:

$$
D^{\mathrm{norm}}_{ij} =
\frac{D^{(1)}_{ij} - \mu_j}{\sqrt{\mu_j}}
$$

where $\mu_j$ is the population mean of $D^{(1)}_{\cdot j}$. Regions are filtered by a variance ratio $\sigma_j^2 / \mu_j$, and the top fraction of high-variance regions are retained for neighbor computation.

**Neighbor-normalized diploid copy number.** GRiD uses nearest-neighbor relationships in the normalized depth space to provide a local reference for copy-number estimation. Following the approach of Mukamel et al. [@mukamel2021], for individual $i$ with read count $r_i$ and depth scale $s_i$, the diploid copy number estimate is:

$$
\widehat{\mathrm{dipCN}}_i =
\frac{r_i/s_i}
{\frac{1}{|N_i|}\sum_{j \in N_i} r_j/s_j}
$$

where $N_i$ is the set of nearest neighbors identified from the normalized depth space using the nearest-neighbor implementation provided by scikit-learn [@pedregosa2011].

**Haplotype inference.** Haplotype-specific copy-number inference follows the approach of Hujoel et al. [@hujoel2026]. This step uses identity-by-descent (IBD) relationships between samples to partition the diploid copy-number estimate into haplotype-specific estimates. For each individual, the two haplotype estimates are iteratively updated according to the copy numbers of IBD-matched haplotypes:

$$
h_{i,k}^{(t+1)} =
\widehat{\mathrm{dipCN}}_i
\cdot
\frac{\bar{h}_{N_{i,k}}^{(t)}}
{\bar{h}_{N_{i,1}}^{(t)} + \bar{h}_{N_{i,2}}^{(t)}}
$$

where $\bar{h}_{N_{i,k}}^{(t)}$ is the mean haplotype copy number of IBD neighbors on haplotype $k$ at iteration $t$. At each iteration, the relative copy numbers of IBD-matched haplotypes determine the allocation of the individual's diploid copy-number estimate between the two haplotypes. The estimates are iteratively updated to incorporate information from IBD-related haplotypes.

# AI Usage Disclosure

The authors used GitHub Copilot and Anthropic Claude (Claude Sonnet 4) during the development of this software and the preparation of this manuscript. GitHub Copilot was used as an interactive coding assistant to suggest code completions, refactoring ideas, and boilerplate implementations. Claude was used to assist with code review, debugging, documentation development, and editorial improvements to the manuscript, including suggestions for clarity, organization, and grammar. AI tools were not used to generate scientific conclusions or interpret results. All AI-assisted code and manuscript edits were reviewed, validated, and, where necessary, modified by the authors before inclusion. The authors take full responsibility for the accuracy, correctness, and content of the software and manuscript.

# Acknowledgements

We are deeply indebted to Drs. Po-Ru Loh and Margaux Hujoel for sharing code and providing additional technical support. We would additionally like to thank Dr. Akshay Avvaru for his help and support. This work was supported by NIH K01: HL157658, R01: HL172887 and NSF: 2022138 grants. The funders had no role in the design, development, analysis, or preparation of this work.

# References
