#!/bin/bash

set -euo pipefail

module load samtools

NUM_NEIGHBORS=200
FOCAL_BP="" # LPA KIV-2 region center, eg 160,626,361
while [[ $# -gt 0 ]]; do
    case $1 in
        --neighbors)
            NUM_NEIGHBORS="$2"; shift 2 ;;
        --focal-bp)
            FOCAL_BP="$2"; shift 2 ;;
        -h|--help)
            echo "Usage: $0 [--neighbors N]"
            exit 0 ;;
        *)
            echo "Unknown argument: $1"
            exit 1 ;;
    esac
done

# Paths and parameters
WORK_DIR="${WORK_DIR:-$(pwd)}"
LOG_DIR="$WORK_DIR/logs"
DATA_DIR="$WORK_DIR/data"
OUTPUT_DIR="$WORK_DIR/output"
mkdir -p "$DATA_DIR" "$LOG_DIR" "$OUTPUT_DIR"
 
THREADS="${SLURM_CPUS_PER_TASK:-$(nproc)}"

timestamp() { date '+%Y-%m-%d %H:%M:%S'; }
log() { echo "[$(timestamp)] $*" | tee -a "$LOG_DIR/IBS.log"; }

download_with_retry() {
    local url="$1"
    local out="$2"
    local attempt
 
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

# Dependency check
MISSING=()
for cmd in qctool wget samtools; do
    command -v "$cmd" &>/dev/null || MISSING+=("$cmd")
done
if [[ ${#MISSING[@]} -gt 0 ]]; then
    log "ERROR: missing required tools: ${MISSING[*]}"
    log "  qctool:   https://www.well.ox.ac.uk/~gav/qctool_v2/"
    log "  samtools: conda install -c bioconda samtools"
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

# 1000 Genomes public URLs
PHASED_VCF_URL="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr6.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
GENETIC_MAP_URL="https://alkesgroup.broadinstitute.org/Eagle/downloads/tables/genetic_map_hg38_withX.txt.gz"
REGIONS_FILE_URL="https://raw.githubusercontent.com/caterer-z-t/GRiD/main/files/734_possible_coding_vntr_regions.IBD2R_gt_0.25.uniq.txt"

# Parameters
# Fetch static inputs
if [[ ! -f "$DATA_DIR/regions.txt" ]]; then
    log "Downloading VNTR regions file..."
    download_with_retry \
        "$REGIONS_FILE_URL" \
        "$DATA_DIR/regions.txt"
fi
 
read -r CHR START END < <(awk '$7=="LPA" {print $1, $2, $3; exit}' "$DATA_DIR/regions.txt") || true
if [[ -z "${CHR:-}" || -z "${START:-}" || -z "${END:-}" ]]; then
    log "ERROR: Could not parse LPA coordinates from regions.txt"
    exit 1
fi
REGION="chr${CHR}:${START}-${END}"

PHASED_VCF="$DATA_DIR/1kGP_chr6_phased.vcf.gz"
LPA_VCF="$DATA_DIR/1kGP_chr6_LPA_phased.vcf.gz"
BGEN_FILE="$DATA_DIR/1kGP_chr6_LPA_phased.bgen"
SAMPLE_FILE="$DATA_DIR/1kGP_chr6_LPA_phased.sample"
GENETIC_MAP="$DATA_DIR/genetic_map_hg38_chr6.txt"
OUTPUT_FILE="$DATA_DIR/ibs_neighbors_chr6.tsv.gz"

# Download 1000G phased VCF for chr6
if [[ ! -f "$LPA_VCF" ]]; then
    log "Streaming LPA region from 1000G phased VCF..."

    # Fetch index alongside to enable remote tabix streaming
    download_with_retry \
        "${PHASED_VCF_URL}.tbi" \
        "${PHASED_VCF}.tbi"

    # Pull only the LPA window
    tabix -h "$PHASED_VCF_URL" "$REGION" \
        | bgzip > "$LPA_VCF"
    tabix -p vcf "$LPA_VCF"

    rm -f "${PHASED_VCF}.tbi"

    log "LPA phased VCF ready: $LPA_VCF"
else
    log "LPA phased VCF already present — skipping."
fi

# Phase 2 — Convert phased VCF → BGEN v1.2
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

# Download Eagle genetic map for hg38
if [[ ! -f "$GENETIC_MAP" ]]; then
    log "Downloading Eagle hg38 genetic map..."
    download_with_retry \
        "$GENETIC_MAP_URL" \
        "$GENETIC_MAP"
    log "Genetic map ready: $GENETIC_MAP"
else
    log "Genetic map already present — skipping."
fi

# Phase 4 — Run computeIBSpbwt
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
