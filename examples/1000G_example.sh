#!/bin/bash
# =============================================================================
# GRiD — 1000 Genomes Project Example
# LPA KIV-2 diploid copy number estimation (hg38, chr6)
#
# Uses publicly available 30x WGS CRAMs from the 1000G high-coverage dataset.
# CRAMs are streamed remotely — only the LPA region is pulled per sample,
# so this does not require downloading 15 GB files.
#
# Usage:
#   bash 1000G_example.sh                  # run all samples (2504)
#   bash 1000G_example.sh --pop EUR        # filter to a superpopulation
#   bash 1000G_example.sh --n 50           # limit to first N samples
#   bash 1000G_example.sh --pop EUR --n 50
#
# Runtime: ~30 min for 100 samples on 8 threads (streaming-bound)
# =============================================================================

set -euo pipefail

# -----------------------------------------------------------------------
# Parse arguments
# -----------------------------------------------------------------------
POP_FILTER=""
N_SAMPLES=0   # 0 = all

while [[ $# -gt 0 ]]; do
    case $1 in
        --pop) POP_FILTER="$2"; shift 2 ;;
        --n)   N_SAMPLES="$2";  shift 2 ;;
        *) echo "Unknown argument: $1"; exit 1 ;;
    esac
done

# -----------------------------------------------------------------------
# Dependency check
# -----------------------------------------------------------------------
MISSING=()
for cmd in samtools mosdepth conda awk wget; do
    command -v "$cmd" &>/dev/null || MISSING+=("$cmd")
done
if [[ ${#MISSING[@]} -gt 0 ]]; then
    echo "ERROR: missing required tools: ${MISSING[*]}"
    echo "Install with: conda install -c bioconda ${MISSING[*]}"
    exit 1
fi

# -----------------------------------------------------------------------
# Paths and parameters
# -----------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
GRID_ROOT="$(dirname "$SCRIPT_DIR")"                      # repo root

WORK_DIR="${WORK_DIR:-$(pwd)/1000G_grid_run}"             # override with env var
CRAM_DIR="$WORK_DIR/crams"
LOG_DIR="$WORK_DIR/logs"
DATA_DIR="$WORK_DIR/data"
OUTPUT_DIR="$WORK_DIR/output"
MOSDEPTH_WORK="$WORK_DIR/mosdepth_work"

THREADS="${SLURM_CPUS_PER_TASK:-$(nproc)}"                # auto-detect cores
STREAM_JOBS=$(( THREADS > 4 ? 4 : THREADS ))              # parallel CRAM streams

# -----------------------------------------------------------------------
# 1000 Genomes public URLs (no login required)
# -----------------------------------------------------------------------
BASE_URL="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp"
CRAM_BASE="${BASE_URL}/data_collections/1000G_2504_high_coverage/working/20190425_NYGC_GATK"
PANEL_URL="${BASE_URL}/release/20130502/integrated_call_samples_v3.20200731.ALL.panel"
REF_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"

# VNTR region of interest — read from bundled regions file
REGIONS_FILE="$GRID_ROOT/files/734_possible_coding_vntr_regions.IBD2R_gt_0.25.uniq.txt"
read -r CHR START END < <(awk '$7=="LPA" {print $1, $2, $3; exit}' "$REGIONS_FILE")
REGION="chr${CHR}:${START}-${END}"

REPEAT_MASK="$GRID_ROOT/files/repeat_mask_list.hg38.ucsc_bed"

# -----------------------------------------------------------------------
# Setup
# -----------------------------------------------------------------------
mkdir -p "$CRAM_DIR" "$LOG_DIR" "$DATA_DIR" "$OUTPUT_DIR" "$MOSDEPTH_WORK"

timestamp() { date '+%Y-%m-%d %H:%M:%S'; }
log() { echo "[$(timestamp)] $*" | tee -a "$LOG_DIR/pipeline.log"; }

log "================================================="
log "  GRiD — 1000 Genomes Example"
log "  Locus:   LPA KIV-2  ($REGION, hg38)"
log "  Threads: $THREADS  |  Parallel streams: $STREAM_JOBS"
log "  Output:  $OUTPUT_DIR"
log "================================================="

# -----------------------------------------------------------------------
# Activate environment
# -----------------------------------------------------------------------
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate grid-env

# -----------------------------------------------------------------------
# Phase 1 — Reference genome
# -----------------------------------------------------------------------
REF_FA="$DATA_DIR/GRCh38_no_alt.fa"

if [[ ! -f "$REF_FA" ]]; then
    log "Downloading GRCh38 reference genome..."
    wget -q -O - "$REF_URL" | gunzip > "$REF_FA"
    samtools faidx "$REF_FA"
    log "Reference genome ready: $REF_FA"
else
    log "Reference genome already present — skipping download."
fi

# -----------------------------------------------------------------------
# Phase 2 — Sample list
# -----------------------------------------------------------------------
PANEL="$DATA_DIR/1000G_panel.txt"
SAMPLES_FILE="$DATA_DIR/samples.txt"

if [[ ! -f "$PANEL" ]]; then
    log "Downloading 1000G panel file..."
    wget -q -O "$PANEL" "$PANEL_URL"
fi

# Columns: sample  pop  super_pop  gender
if [[ -n "$POP_FILTER" ]]; then
    log "Filtering to superpopulation: $POP_FILTER"
    awk -v pop="$POP_FILTER" 'NR>1 && $3==pop {print $1}' "$PANEL" > "$SAMPLES_FILE"
else
    awk 'NR>1 {print $1}' "$PANEL" > "$SAMPLES_FILE"
fi

TOTAL=$(wc -l < "$SAMPLES_FILE")
if [[ "$N_SAMPLES" -gt 0 && "$N_SAMPLES" -lt "$TOTAL" ]]; then
    head -n "$N_SAMPLES" "$SAMPLES_FILE" > "${SAMPLES_FILE}.tmp" && mv "${SAMPLES_FILE}.tmp" "$SAMPLES_FILE"
    log "Using first $N_SAMPLES of $TOTAL available samples."
else
    log "Using all $TOTAL samples."
fi

# -----------------------------------------------------------------------
# Phase 3 — Stream LPA-region CRAMs (parallel, remote)
#
# Only the target region is pulled over the network per sample.
# A 30x whole-genome CRAM is ~15 GB; the LPA region slice is ~2 MB.
# -----------------------------------------------------------------------
log "Streaming LPA region CRAMs from 1000G EBI..."

stream_cram() {
    local sample="$1"
    local url="${CRAM_BASE}/${sample}.hg38.cram"
    local out="${CRAM_DIR}/${sample}.cram"

    if [[ -f "$out" && -f "${out}.crai" ]]; then
        return 0   # already done
    fi

    samtools view -T "$REF_FA" -b "$url" "$REGION" \
        | samtools sort -@ 2 -o "$out" && \
    samtools index "$out"
}

export -f stream_cram
export REF_FA CRAM_BASE CRAM_DIR REGION

# Use GNU parallel if available, otherwise xargs
if command -v parallel &>/dev/null; then
    parallel -j "$STREAM_JOBS" --bar stream_cram :::: "$SAMPLES_FILE" \
        >> "$LOG_DIR/stream_crams.log" 2>&1
else
    xargs -P "$STREAM_JOBS" -I {} bash -c 'stream_cram "$@"' _ {} \
        < "$SAMPLES_FILE" >> "$LOG_DIR/stream_crams.log" 2>&1
fi

log "CRAM streaming complete."

# -----------------------------------------------------------------------
# Phase 4 — Generate config
# -----------------------------------------------------------------------
CONFIG="$WORK_DIR/config.yaml"

cat > "$CONFIG" << YAML
# Auto-generated GRiD config — 1000 Genomes LPA KIV-2 example
samples_file: "$SAMPLES_FILE"
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
  run: False                     # CRAMs already indexed above

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
    min_depth: 20
    max_depth: 100
    top_frac: 0.1
    output_file_prefix: "mosdepth_normalized"
    repeat_mask_file: "$REPEAT_MASK"

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
  run: False                     # requires IBS file — see 1000G_IBS_example.sh
  output_file_prefix: "haploid_genotypes"
  ibs_output: "$DATA_DIR/ibs_neighbors_chr6.tsv.gz"
  min_neighbors: 1
  max_neighbors: 10
  n_iters: 100
YAML

log "Config written to: $CONFIG"

# -----------------------------------------------------------------------
# Phase 5 — Run GRiD
# -----------------------------------------------------------------------
log "Running GRiD pipeline..."
grid wgs "$CONFIG" 2>&1 | tee "$LOG_DIR/grid.log"

log "================================================="
log "  Done. Diploid CN results in: $OUTPUT_DIR"
log "  Next: run 1000G_IBS_example.sh to generate the"
log "  IBS neighbors file, then re-run with"
log "  compute_haploid_genotypes.run: True"
log "================================================="
