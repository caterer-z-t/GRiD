# GRiD — Genomic Repeat inference from Depth

[![Python](https://img.shields.io/badge/python-3.8%2B-blue?logo=python&logoColor=white)](https://www.python.org/)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![pytest](https://img.shields.io/badge/tested%20with-pytest-0A9EDC?logo=pytest&logoColor=white)](https://docs.pytest.org/)
[![Coverage](https://codecov.io/gh/caterer-z-t/GRiD/branch/main/graph/badge.svg)](https://codecov.io/gh/caterer-z-t/GRiD)
[![DOI](https://joss.theoj.org/papers/placeholder/status.svg)](https://joss.theoj.org)
[![GitHub issues](https://img.shields.io/github/issues/caterer-z-t/GRiD)](https://github.com/caterer-z-t/GRiD/issues)
[![GitHub stars](https://img.shields.io/github/stars/caterer-z-t/GRiD?style=social)](https://github.com/caterer-z-t/GRiD/stargazers)
[![Visits](https://visitor-badge.laobi.icu/badge?page_id=caterer-z-t.GRiD)](https://github.com/caterer-z-t/GRiD)

**GRiD** is a Python pipeline for haplotype-resolved VNTR copy number estimation from short-read whole-genome sequencing data. It was developed for and validated on the KIV-2 VNTR in the *LPA* gene, but the core depth-normalization framework is applicable to any large VNTR locus.

The *LPA* KIV-2 repeat is one of the largest and most structurally complex VNTRs in the human genome (spanning up to 200 kb) making standard STR and CNV tools inadequate. GRiD addresses this by combining read-depth normalization, ancestry-matched nearest-neighbor estimation, and IBD-based haplotype inference into a single, configurable pipeline.

---

## Installation

Requires Python ≥ 3.8 and [`mosdepth`](https://github.com/brentp/mosdepth) available on your `PATH`.

```bash
git clone https://github.com/caterer-z-t/GRiD.git
cd GRiD

# Create environment
conda env create -f env.yaml
conda activate grid-env

# Install
pip install -e .
```

---

## Quick Start

GRiD is driven by a single YAML configuration file. Copy and edit the example config:

```bash
cp grid/example_config.yaml my_config.yaml
# edit my_config.yaml with your file paths and parameters
```

Run the full WGS pipeline:

```bash
grid wgs my_config.yaml
```

See [`grid/example_config.yaml`](https://github.com/caterer-z-t/GRiD/tree/main/grid/example_config.yaml) for a fully annotated configuration reference.

---

## Pipeline

GRiD runs seven sequential steps, each controlled by a `run: True/False` flag in the config:

| Step | Module | Description | Locus-specific? |
|------|--------|-------------|-----------------|
| 1 | `ensure_crai` | Verify or create CRAM index files | No |
| 2 | `count_reads` | Count properly paired reads in the target VNTR region | No |
| 3 | `mosdepth` | Compute binned read depth across the genome | No |
| 4 | `normalize_mosdepth` | Normalize coverage within and across individuals | No |
| 5 | `find_neighbors` | Identify read depth-matched nearest neighbors | No |
| 6 | `compute_dipcn` | Estimate diploid copy number via neighbor normalization | No |
| 7 | `hi_inference` | Infer haplotype-specific copy numbers via IBD phasing | No* |

$*$ Step 7 requires a Identity By Descent/State output file from [Hujoel 2026 Nature](https://doi.org/10.1038/s41586-025-09886-z) `.cpp` file; for further information and how to obtain this file, please see [IBS/IBD](https://caterer-z-t.github.io/GRiD/IBS_IBD.html) page of our documentation. 

For full algorithmic details, see the [Algorithms page](https://caterer-z-t.github.io/GRiD/algorithms.html) of our documentation.

---

## Configuration

All parameters are set in a YAML file. Top-level fields define global settings; each step has its own section with a `run` flag and step-specific options.

<div style="padding:10px; border-left:4px solid #FF7800;">
<strong>Note:</strong> Each step is controlled by a <code>run: True/False</code> flag. Setting a step to <code>False</code> skips it entirely — and because the pipeline is sequential, any downstream step that depends on its output will also fail to run correctly. Enable steps in order and ensure all prerequisite steps have been run at least once before skipping them.
</div>


```yaml
samples_file: "path/to/sample_ids.txt"
directory_loc: "path/to/crams/"
reference_genome: "path/to/hg38.fa"
output_dir: "path/to/output/"
threads: 8
file_type: "cram"
chrom: "chr6"
start_bp: 160605062
end_bp: 160647661

mosdepth:
  run: True
  bin_size: 1000
  ...

compute_haploid_genotypes:
  run: True
  ibs_output: "path/to/ibs_neighbors_chr6.tsv.gz"
  min_neighbors: 1
  max_neighbors: 10
  n_iters: 100
```

See [`grid/example_config.yaml`](./grid/example_config.yaml) for all available options and their defaults.

For a full breakdown of every parameter, how to obtain required input files, and configuration tips see the [Getting Started page](https://caterer-z-t.github.io/GRiD/getting_started.html).

---

## HPC Usage

For large cohorts on SLURM-based clusters, see [`examples/SOL.sh`](./examples/1000G_example.sh) for an example submission script used on the 1000 Genomes on University of North Carolina's Super Computing Cluster Longleaf.

---

## Citation

If you use GRiD, please cite:

```bibtex
@software{caterer2026,
    author  = {Caterer, Zachary, Lin Meng, An Qiang, Lange Ethan, Gignoux Christopher, Graff Mariaelisa, Avery Christy, Stanislawski Maggie},
    title   = {GRiD: Genomic Repeat Inference from Depth},
    year    = {2026},
    url     = {https://github.com/caterer-z-t/GRiD},
    version = {1.0.0}
}
```

GRiD builds on methods from:
- Mukamel et al. (2021) [doi:10.1126/science.abg8635](https://doi.org/10.1126/science.abg8635)
- Hujoel et al. (2026) [doi:10.1038/s41586-025-09886-z](https://doi.org/10.1038/s41586-025-09886-z)

---

## License

Released under the [MIT License](./LICENSE).

## Support

- [Issues](https://github.com/caterer-z-t/GRiD/issues)
- [Discussions](https://github.com/caterer-z-t/GRiD/discussions)
- [Contributing](./CONTRIBUTING.md)
- [Code of Conduct](./CODE_OF_CONDUCT.md)
