#!/bin/bash
#SBATCH --job-name=grid_run
#SBATCH --output=slurm/%x.out
#SBATCH --error=slurm/%x.err
#SBATCH --time=4-00:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G

# Stage 3 of 3 — generate config.yaml and run the GRiD pipeline.
# Submit after 01_prep_data.sh (and, if you want haploid CN, 02_run_ibs.sh)
# have both completed.
#
# Usage: sbatch 03_generate_config_and_run.sh   (from the same working
# directory as stages 1 and 2)

set -euo pipefail

module load anaconda mosdepth samtools

WORK_DIR="$(pwd)"
CRAM_DIR="$WORK_DIR/crams"
LOG_DIR="$WORK_DIR/logs"
DATA_DIR="$WORK_DIR/data"
OUTPUT_DIR="$WORK_DIR/output"
MOSDEPTH_WORK="$WORK_DIR/mosdepth_work"
mkdir -p "$OUTPUT_DIR" "$MOSDEPTH_WORK" "$LOG_DIR" "$WORK_DIR/slurm"

THREADS="${SLURM_CPUS_PER_TASK:-$(nproc)}"

timestamp() { date '+%Y-%m-%d %H:%M:%S'; }
log() { echo "[$(timestamp)] $*" | tee -a "$LOG_DIR/pipeline.log"; }

# --- Load locus coordinates written by 01_prep_data.sh ---
VARS_FILE="$DATA_DIR/pipeline_vars.sh"
if [[ ! -f "$VARS_FILE" ]]; then
    echo "ERROR: $VARS_FILE not found — run 01_prep_data.sh first."
    exit 1
fi
source "$VARS_FILE"

REF_FA="$DATA_DIR/GRCh38_no_alt.fa"
REPEAT_MASK="$DATA_DIR/repeat_mask.bed"
GRID_SAMPLES_FILE="$DATA_DIR/grid_samples.txt"

for f in "$REF_FA" "$GRID_SAMPLES_FILE"; do
    if [[ ! -f "$f" ]]; then
        echo "ERROR: expected input $f not found — run 01_prep_data.sh first."
        exit 1
    fi
done

# --- Activate GRiD environment ---
if ! conda env list | grep -qE '^\s*grid-test\s'; then
    echo "ERROR: conda environment 'grid-test' not found."
    echo "Create it first, e.g.: conda create -n grid-test python=3.11 -y"
    exit 1
fi
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate grid-test

if ! command -v grid &>/dev/null; then
    echo "ERROR: 'grid' command not found after activating grid-test."
    echo "Check that 'pip install -e .' completed successfully in that env."
    exit 1
fi

# --- Check for IBS output from stage 2 (optional — diploid-only if absent) ---
IBS_OUTPUT="$DATA_DIR/ibs_neighbors_chr${CHR}.tsv.gz"
HAPLOID_RUN="False"
if [[ -f "$IBS_OUTPUT" ]]; then
    log "IBS neighbors file found: $IBS_OUTPUT"
    log "Haploid CN estimation will run."
    HAPLOID_RUN="True"
else
    log "IBS neighbors file not found at $IBS_OUTPUT."
    log "Run 02_run_ibs.sh first if you want haploid CN estimation."
    log "Proceeding with diploid-only for this run."
fi

# Empty mask: the downloaded repeat_mask.bed spans nearly the entire LPA
# window (KIV-2 is itself a massive tandem repeat), which zeroes out every
# normalization bin if applied here. See prior discussion — using an empty
# mask for this locus is the correct choice, not a workaround being skipped.
EMPTY_MASK="$DATA_DIR/repeat_mask_empty.bed"
touch "$EMPTY_MASK"

# --- Generate config ---
CONFIG="$WORK_DIR/config.yaml"

cat > "$CONFIG" << YAML
# Auto-generated GRiD config — 1000 Genomes LPA KIV-2 example
samples_file: "$GRID_SAMPLES_FILE"
directory_loc: "$CRAM_DIR"
reference_genome: "$REF_FA"
output_dir: "$OUTPUT_DIR"
threads: $THREADS
file_type: "cram"
chrom: "chr${CHR}"
start_bp: $START
end_bp: $END
output_file_type: "tsv"

index:
  run: True
  output_file_prefix: "index_file_results"

count_reads:
  run: True
  output_file_prefix: "read_counts"
  min_mapq: 1
  flags:
    - 83    # proper pair, read reverse strand
    - 147   # proper pair, mate reverse strand
    - 81    # read reverse strand
    - 145   # mate reverse strand

mosdepth:
  run: True
  output_file_prefix: "mosdepth_results"
  bin_size: 1000
  mode: "fast"
  region_name: "LPA"
  work_dir: "$MOSDEPTH_WORK"
  remove_intermediate: True

  normalize:
    run: True
    min_depth: 1
    max_depth: 30
    top_frac: 0.9
    output_file_prefix: "mosdepth_normalized"
    repeat_mask_file: "$EMPTY_MASK"

  neighbors:
    run: True
    output_file_prefix: "neighbors"
    num_neighbors: 5
    zmax: 2.0
    sigma2_max: 1000

compute_diploid_genotypes:
  run: True
  output_file_prefix: "diploid_genotypes"

compute_haploid_genotypes:
  run: $HAPLOID_RUN
  output_file_prefix: "haploid_genotypes"
  ibs_output: "$IBS_OUTPUT"
  min_neighbors: 1
  max_neighbors: 10
  n_iters: 100
  method: "ibs"
YAML

log "Config written to: $CONFIG"
log "  compute_haploid_genotypes.run: $HAPLOID_RUN"

# --- Run GRiD ---
log "Running GRiD pipeline..."
grid wgs "$CONFIG"

log "Done. Results in: $OUTPUT_DIR"