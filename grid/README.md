# GRiD Core Module

This directory contains the core implementation of the GRiD pipeline for LPA KIV-2 copy number variant estimation.

## Structure

```
grid/
├── cli.py              # Command-line interface with all pipeline commands
├── config.py           # Configuration management (in development)
├── config.yaml         # Example configuration file
├── pipeline.py         # Unified pipeline orchestration (in development)
├── utils/              # Core computational modules and utilities
└── README.md           # This file
```

## Key Components

### `cli.py`
The main entry point for all GRiD commands. Provides a rich CLI interface with commands for:
- CRAM file indexing and subsetting
- Read counting in VNTR regions
- Coverage analysis with mosdepth
- Neighbor finding and normalization
- Reference extraction
- LPA realignment
- Copy number estimation

Run `grid --help` to see all available commands.

### `config.yaml`
Example configuration file for setting up genomic coordinates, file paths, and pipeline parameters.

### `utils/`
Contains all computational logic organized into modular Python scripts and subdirectories. See `utils/README.md` for detailed module information.

## Usage

All commands are accessed through the `grid` CLI:

```bash
# Get help
grid --help
grid <command> --help

# Example: Count reads
grid count-reads --cram-dir data/CRAMs \
                 --output-file output/counts.tsv \
                 --ref-fasta refs/hs37d5.fa \
                 --chrom chr6 \
                 --start 160000000 \
                 --end 160100000 \
                 --config config.yaml
```

## Future Development

- [ ] Implement `pipeline.py` for unified `grid run` command
- [ ] Enhance `config.py` for advanced configuration management
- [ ] Add Python API for programmatic access