#!/bin/bash
#SBATCH --job-name=grid_ibs
#SBATCH --output=slurm/%x.out
#SBATCH --error=slurm/%x.err
#SBATCH --time=2-00:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=100G

# Stage 2 of 3 — compute IBS neighbors (needed for haploid CN estimation).
# Submit after 01_prep_data.sh completes; sources data/pipeline_vars.sh for
# CHR/FOCAL_BP rather than re-deriving them, so this can't disagree with
# stage 1's locus coordinates.
#
# Usage: sbatch 02_run_ibs.sh   (from the same working directory as stage 1)

set -euo pipefail

module load samtools

NUM_NEIGHBORS=200
while [[ $# -gt 0 ]]; do
    case $1 in
        --num-neighbors) NUM_NEIGHBORS="$2"; shift 2 ;;
        -h|--help)
            echo "Usage: $0 [--num-neighbors N]"
            exit 0 ;;
        *)
            echo "Unknown argument: $1"
            exit 1 ;;
    esac
done

WORK_DIR="${WORK_DIR:-$(pwd)}"
LOG_DIR="$WORK_DIR/logs"
DATA_DIR="$WORK_DIR/data"
OUTPUT_DIR="$WORK_DIR/output"
mkdir -p "$DATA_DIR" "$LOG_DIR" "$OUTPUT_DIR" "$WORK_DIR/slurm"

THREADS="${SLURM_CPUS_PER_TASK:-$(nproc)}"

timestamp() { date '+%Y-%m-%d %H:%M:%S'; }
log() { echo "[$(timestamp)] $*" | tee -a "$LOG_DIR/IBS.log"; }

download_with_retry() {
    local url="$1" out="$2" attempt
    for attempt in 1 2 3 4 5; do
        if wget -q --tries=1 --timeout=60 -O "$out" "$url"; then
            return 0
        fi
        rm -f "$out"
        sleep $(( attempt * 3 + RANDOM % 3 ))
    done
    echo "ERROR: failed to download $url after 5 attempts" >&2
    return 1
}

# --- Load locus coordinates written by 01_prep_data.sh ---
VARS_FILE="$DATA_DIR/pipeline_vars.sh"
if [[ ! -f "$VARS_FILE" ]]; then
    echo "ERROR: $VARS_FILE not found — run 01_prep_data.sh first."
    exit 1
fi
source "$VARS_FILE"

if ! command -v "$WORK_DIR/qctool/qctool" &> /dev/null; then
    log "Installing qctool..."
    download_with_retry \
        "https://www.well.ox.ac.uk/~gav/resources/qctool_v2.2.0-CentOS_Linux7.8.2003-x86_64.tgz" \
        "$WORK_DIR/qctool.tgz"
    tar -xzf "$WORK_DIR/qctool.tgz" -C "$WORK_DIR"
    mv "$WORK_DIR/qctool_v2.2.0-CentOS Linux7.8.2003-x86_64" "$WORK_DIR/qctool"
    rm -f "$WORK_DIR/qctool.tgz"
fi

MISSING=()
for cmd in "$WORK_DIR/qctool/qctool" wget samtools; do
    command -v "$cmd" &>/dev/null || MISSING+=("$cmd")
done
if [[ ${#MISSING[@]} -gt 0 ]]; then
    log "ERROR: missing required tools: ${MISSING[*]}"
    exit 1
fi

COMPUTE_IBS="$WORK_DIR/computeIBSpbwt"
if [[ ! -f "$COMPUTE_IBS" ]]; then
    log "Downloading computeIBSpbwt binary..."
    download_with_retry \
        "https://raw.githubusercontent.com/mhujoel/STRs/main/cpp_files/computeIBSpbwt" \
        "$COMPUTE_IBS"
    chmod +x "$COMPUTE_IBS"
fi

PHASED_VCF_URL="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr${CHR}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
GENETIC_MAP_URL="https://alkesgroup.broadinstitute.org/Eagle/downloads/tables/genetic_map_hg38_withX.txt.gz"

LPA_VCF="$DATA_DIR/1kGP_chr${CHR}_LPA_phased.vcf.gz"
BGEN_FILE="$DATA_DIR/1kGP_chr${CHR}_LPA_phased.bgen"
SAMPLE_FILE="$DATA_DIR/1kGP_chr${CHR}_LPA_phased.sample"
GENETIC_MAP="$DATA_DIR/genetic_map_hg38_chr${CHR}.txt"
OUTPUT_FILE="$DATA_DIR/ibs_neighbors_chr${CHR}.tsv.gz"

# --- Stream LPA region from the 1000G phased VCF ---
if [[ ! -f "$LPA_VCF" ]]; then
    log "Streaming LPA region from 1000G phased VCF..."
    tabix -h "$PHASED_VCF_URL" "$REGION" | bgzip > "$LPA_VCF"
    tabix -p vcf "$LPA_VCF"
    find "$WORK_DIR" -maxdepth 1 -name "*.tbi" ! -path "${LPA_VCF}.tbi" -delete
    log "LPA phased VCF ready: $LPA_VCF"
else
    log "LPA phased VCF already present — skipping."
fi

# --- Convert phased VCF -> BGEN v1.2 ---
if [[ ! -f "$BGEN_FILE" ]]; then
    log "Converting phased VCF to BGEN v1.2..."
    "$WORK_DIR/qctool/qctool" \
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

# --- Download + decompress + filter Eagle genetic map ---
# The map's chrom column is bare (e.g. "6"), not "chr6" — filtering on
# "chr${CHR}" here would silently match nothing and produce a header-only
# file, so this must stay as the bare chromosome number.
if [[ ! -f "$GENETIC_MAP" ]]; then
    log "Downloading Eagle hg38 genetic map..."
    wget -q -O - "$GENETIC_MAP_URL" \
        | gunzip \
        | awk -v c="${CHR}" 'NR==1 || $1==c' \
        > "$GENETIC_MAP"
    log "Genetic map ready: $GENETIC_MAP"
else
    log "Genetic map already present — skipping."
fi

# --- Run computeIBSpbwt ---
log "Computing IBS neighbors..."
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

log "IBS neighbors written to: $OUTPUT_FILE"
log "Next: submit 03_generate_config_and_run.sh"