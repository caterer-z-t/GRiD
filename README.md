# GRiD - Genomic Repeat inference from Depth

## LPA KIV-2 Copy Number Variant Estimation Pipeline

**GRiD** is a modular, extensible pipeline for estimating copy number variants (CNVs) in the *LPA* gene’s **KIV-2 VNTR** region using short-read sequencing data.

Although currently specialized for *LPA*, GRiD is designed with modularity in mind — future extensions to other repeat loci are encouraged.

Based on the methodology from:
- [Mukamel et al. (2021)](https://doi.org/10.1126/science.abg8289) - Original pipeline
- [Hao (2024)](https://github.com/alexliyihao/vntrwrap) - C++/shell implementation


## Installation

It is strongly recommended to use a dedicated environment (e.g., `conda`, `micromamba`, or `venv`). 

### From source
```bash
git clone https://github.com/caterer-z-t/GRiD.git
cd GRiD

# Create environment (recommended)
conda env create -n grid-env -f env.yaml
conda activate grid-env

# Install in editable mode (for now, we will incoorperate this into pypi + conda)
pip install -e .
```

## Usage

GRiD provides a modular **command-line interface (CLI)** for reproducible and scalable execution on Linux or HPC environments (e.g., UNC Longleaf).

### Basic Example
```bash
grid count-reads --cram-dir data/CRAMs \
                 --output-file output/counts.tsv \
                 --ref-fasta refs/hs37d5.fa \
                 --chrom chr6 \
                 --start 160000000 \
                 --end 160100000 \
                 --config config.yaml
```

### Basic Pipeline Run

For large CRAM datasets, the CLI is the recommended interface.
See [examples/SOL.sh](https://caterer-z-t.github.io/GRiD/examples/SOL_sbatch.html) for an example SLURM submission script.

    
> **Note:** 
The included `SOL.sh` script demonstrates use on the **Study of Latinas (SOL)** dataset, subset via a [Terra](https://app.terra.bio/) workflow and executed on UNC Longleaf (Red Hat 9.6). 
See [Longleaf](https://help.rc.unc.edu/longleaf-cluster/) documentation for environment-specific details.

### Configuration

GRiD uses YAML configuration files for flexible setup.

Customize paths, genomic coordinates, and environment parameters in your `config.yaml` before execution.

## Pipeline Steps

| Step                      | Description                                                 |
| ------------------------- | ----------------------------------------------------------- |
| 1. **GCloud Copy**        | *(Deprecated)* Copy files from Google Cloud Storage.        |
| 2. **Count Reads**        | Count reads within VNTR regions for each CRAM.              |
| 3. **Mosdepth Coverage**  | Compute read coverage using `mosdepth`.                     |
| 4. **Normalize Coverage** | Normalize per-sample coverage values.                       |
| 5. **Find Neighbors**     | Identify similar samples for reference-based normalization. |
| 6. **Extract Reference**  | Extract VNTR reference sequences.                           |
| 7. **Realign Reads**      | Realign reads to the extracted VNTR reference.              |
| 8. **Compute Diploid CN** | Estimate diploid copy numbers.                              |
| 9. **Estimate KIV CN**    | Estimate KIV-2 copy number per individual.                  |

Each step corresponds to a CLI command (e.g., `grid mosdepth`, `grid normalize-mosdepth`, etc.).

Run `grid --help` for a full command list.

Additionally, each of the corresponding commands also have inline documentation, run `grid mosdepth --help`, `grid normalize-mosdepth --help` for stepwise parameters. 

## Pipeline Orchestration
The `grid/utils/` directory contains all computational logic and helper scripts.

| Module                       | Description                                    |
| ---------------------------- | ---------------------------------------------- |
| `ensure_crai.py`             | Ensures CRAM index (.crai) files exist         |
| `count_reads.py`             | Counts paired reads in VNTR region             |
| `run_mosdepth.py`            | Computes coverage using `mosdepth`             |
| `run_normalized_mosdepth.py` | Normalizes coverage across samples             |
| `find_neighbors.py`          | Identifies nearest neighbors for normalization |
| `extract_reference.py`       | Extracts VNTR reference sequences              |
| `align_lpa.py`               | Realigns reads to *LPA* reference              |
| `compute_dipcn.py`           | Computes diploid copy number estimates         |
| `estimate_kiv.py`            | Derives KIV-2 CN estimates                     |

Each utility can be imported directly in Python, e.g.:

```python
from grid.utils.run_mosdepth import run_mosdepth

run_mosdepth(
    cram_dir="data/crams",
    output_file="out/mosdepth.tsv",
    ref_fasta="refs/hs37d5.fa",
    chrom="chr6",
    start=160000000,
    end=160100000
)
```

## Development

### Code Style
We follow standard Python formatting and type-checking conventions:
```bash
black grid/
isort grid/
flake8 grid/
mypy grid/
```

### Build Documentation

```bash
conda activate grid-env
pip install -r docs/requirements.txt
make -C docs html
```

Documentation is built with Sphinx and published via ReadTheDocs:
    
> [caterer-z-t.github.io/GRiD](https://caterer-z-t.github.io/GRiD)
## Citation

If you use this pipeline, please cite:

- Mukamel et al. (2021) [Original methodology](https://doi.org/10.1126/science.abg8289)
- Hao et al. (2024) [Implementation](https://github.com/alexliyihao/vntrwrap)
- This pipeline [Your paper when published](https://caterer-z-t.github.io/GRiD/CITATION.cff)

## License

Released under [MIT License](https://caterer-z-t.github.io/GRiD/LICENSE)

## Support

- [Documentation](https://caterer-z-t.github.io/GRiD)
- [Issues](https://github.com/caterer-z-t/GRiD/issues)
- [Discussions](https://github.com/caterer-z-t/GRiD/discussions)
- [Contributing](https://caterer-z-t.github.io/GRiD/contributing.html)
- [Code of Conduct](https://caterer-z-t.github.io/GRiD/code_of_conduct.html)

# To-Do's
- [] Add/Implement `pipeline.py` and `config.py` for unified `grid run [args]` entrypoint
- [] Finish README Documentation
- [] Add test suite 
- [] CI/CD workflows

