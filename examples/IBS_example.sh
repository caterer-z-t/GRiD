#!/bin/bash
# =============================================================================
# GRiD — IBS/IBD Neighbor Computation
# LPA KIV-2 (hg38, chr6) using 1000 Genomes phased panel
#
# Pre-requisite for GRiD Step 7 (haplotype inference).
# Runs computeIBSpbwt from Hujoel et al. (2026) Nature.
# Source: https://github.com/mhujoel/STRs/blob/main/cpp_files/computeIBSpbwt.cpp
#
# This script:
#   1. Downloads the 1000G phased VCF for chr6 (no login required)
#   2. Extracts the LPA region and converts to BGEN v1.2 via qctool
#   3. Downloads the Eagle hg38 genetic map
#   4. Runs computeIBSpbwt to produce the IBS neighbors file
#
# Usage:
#   bash   IBS_example.sh          # local run
#   sbatch IBS_example.sh          # SLURM submission
#
# Output: ibs_neighbors_chr6.tsv.gz
#   Pass to GRiD via compute_haploid_genotypes.ibs_output in your config.
# =============================================================================
#SBATCH --job-name=grid_IBS
#SBATCH --output=logs/IBS_%j.out
#SBATCH --error=logs/IBS_%j.err
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G

set -euo pipefail

# -----------------------------------------------------------------------
# Dependency check
# -----------------------------------------------------------------------
MISSING=()
for cmd in qctool wget samtools; do
    command -v "$cmd" &>/dev/null || MISSING+=("$cmd")
done
if [[ ${#MISSING[@]} -gt 0 ]]; then
    echo "ERROR: missing required tools: ${MISSING[*]}"
    echo "  qctool:   https://www.well.ox.ac.uk/~gav/qctool_v2/"
    echo "  samtools: conda install -c bioconda samtools"
    exit 1
fi

# -----------------------------------------------------------------------
# Paths — set WORK_DIR or override individual variables
# -----------------------------------------------------------------------
WORK_DIR="${WORK_DIR:-$(pwd)/1000G_grid_run}"     # match 1000G_example.sh

DATA_DIR="$WORK_DIR/data"
LOG_DIR="$WORK_DIR/logs"
mkdir -p "$DATA_DIR" "$LOG_DIR"

COMPUTE_IBS="${COMPUTE_IBS:-$(pwd)/computeIBSpbwt}"   # compiled binary path

# -----------------------------------------------------------------------
# 1000 Genomes public URLs
# -----------------------------------------------------------------------
PHASED_VCF_URL="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr6.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
GENETIC_MAP_URL="https://alkesgroup.broadinstitute.org/Eagle/downloads/tables/genetic_map_hg38_withX.txt.gz"

# -----------------------------------------------------------------------
# Parameters
# -----------------------------------------------------------------------
CHR=6
FOCAL_BP=160626361       # center of LPA KIV-2 VNTR region (hg38)
LPA_REGION="chr6:160605062-160647661"
NUM_NEIGHBORS=200
THREADS="${SLURM_CPUS_PER_TASK:-$(nproc)}"

PHASED_VCF="$DATA_DIR/1kGP_chr6_phased.vcf.gz"
LPA_VCF="$DATA_DIR/1kGP_chr6_LPA_phased.vcf.gz"
BGEN_FILE="$DATA_DIR/1kGP_chr6_LPA_phased.bgen"
SAMPLE_FILE="$DATA_DIR/1kGP_chr6_LPA_phased.sample"
GENETIC_MAP="$DATA_DIR/genetic_map_hg38_chr6.txt"
OUTPUT_FILE="$DATA_DIR/ibs_neighbors_chr6.tsv.gz"

timestamp() { date '+%Y-%m-%d %H:%M:%S'; }
log() { echo "[$(timestamp)] $*" | tee -a "$LOG_DIR/IBS.log"; }

log "================================================="
log "  IBS/IBD Neighbor Computation"
log "  Locus: LPA KIV-2 ($LPA_REGION, hg38)"
log "  Threads: $THREADS"
log "================================================="

# Verify computeIBSpbwt binary
if [[ ! -x "$COMPUTE_IBS" ]]; then
    echo "ERROR: computeIBSpbwt binary not found or not executable at: $COMPUTE_IBS"
    echo "  Download source: https://github.com/mhujoel/STRs/blob/main/cpp_files/computeIBSpbwt.cpp"
    echo "  Then compile and set: export COMPUTE_IBS=/path/to/computeIBSpbwt"
    exit 1
fi

# -----------------------------------------------------------------------
# Phase 1 — Download 1000G phased VCF for chr6
#
# The full chr6 phased panel is ~5 GB. We stream only the LPA window
# into a local subset VCF using tabix (no full download required).
# -----------------------------------------------------------------------
if [[ ! -f "$LPA_VCF" ]]; then
    log "Streaming LPA region from 1000G phased VCF..."

    # Fetch index alongside to enable remote tabix streaming
    wget -q -c -O "${PHASED_VCF}.tbi" "${PHASED_VCF_URL}.tbi"

    # Pull only the LPA window — ~1 MB vs 5 GB for full chr6
    tabix -h "$PHASED_VCF_URL" "$LPA_REGION" \
        | bgzip > "$LPA_VCF"
    tabix -p vcf "$LPA_VCF"

    log "LPA phased VCF ready: $LPA_VCF"
else
    log "LPA phased VCF already present — skipping."
fi

# -----------------------------------------------------------------------
# Phase 2 — Convert phased VCF → BGEN v1.2
#
# qctool outputs BGEN v1.2, layout 2, 16-bit, phased — the format
# required by computeIBSpbwt.
# -----------------------------------------------------------------------
if [[ ! -f "$BGEN_FILE" ]]; then
    log "Converting phased VCF to BGEN v1.2..."
    qctool \
        -g "$LPA_VCF" \
        -og "$BGEN_FILE" \
        -os "$SAMPLE_FILE" \
        -ofiletype bgen_v1.2 \
        -bgen-bits 16 \
        -bgen-permitted-input-rounding-error 0 \
        >> "$LOG_DIR/qctool.log" 2>&1
    log "BGEN file ready: $BGEN_FILE"
else
    log "BGEN file already present — skipping."
fi

# -----------------------------------------------------------------------
# Phase 3 — Download Eagle genetic map for hg38
# -----------------------------------------------------------------------
if [[ ! -f "$GENETIC_MAP" ]]; then
    log "Downloading Eagle hg38 genetic map..."
    wget -q -O - "$GENETIC_MAP_URL" \
        | gunzip \
        | awk 'NR==1 || $1=="chr6"' \
        > "$GENETIC_MAP"
    log "Genetic map ready: $GENETIC_MAP"
else
    log "Genetic map already present — skipping."
fi

# -----------------------------------------------------------------------
# Phase 4 — Run computeIBSpbwt
# -----------------------------------------------------------------------
log "Computing IBS/IBD neighbors..."
log "  CHR=$CHR  FOCAL_BP=$FOCAL_BP  NEIGHBORS=$NUM_NEIGHBORS  THREADS=$THREADS"

"$COMPUTE_IBS" \
    "$CHR"          \
    "$FOCAL_BP"     \
    "$BGEN_FILE"    \
    "$SAMPLE_FILE"  \
    "$GENETIC_MAP"  \
    "$NUM_NEIGHBORS"\
    "$THREADS"      \
    "$OUTPUT_FILE"

log "================================================="
log "  Done. IBS neighbors written to:"
log "    $OUTPUT_FILE"
log ""
log "  Add to your GRiD config:"
log "    compute_haploid_genotypes:"
log "      run: True"
log "      ibs_output: \"$OUTPUT_FILE\""
log "================================================="
