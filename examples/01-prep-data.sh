#!/bin/bash
#SBATCH --job-name=grid_prep
#SBATCH --output=slurm/%x.out
#SBATCH --error=slurm/%x.err
#SBATCH --time=4-00:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G

# Stage 1 of 3 — download reference/panel/regions/repeat-mask, stream LPA
# CRAMs for all 1000G samples, and build the GRiD-facing samples file.
#
# Submit this first. Writes data/pipeline_vars.sh, which stages 2 and 3
# source to pick up the LPA locus coordinates without re-deriving them.
#
# Usage: sbatch 01_prep_data.sh   (run from the shared working directory —
# stages 2 and 3 must be submitted from this same directory)

set -euo pipefail

module load samtools

WORK_DIR="$(pwd)"
CRAM_DIR="$WORK_DIR/crams"
LOG_DIR="$WORK_DIR/logs"
DATA_DIR="$WORK_DIR/data"

mkdir -p "$CRAM_DIR" "$LOG_DIR" "$DATA_DIR" "$WORK_DIR/slurm"

THREADS="${SLURM_CPUS_PER_TASK:-$(nproc)}"
# Cap concurrent EBI connections — higher concurrency has caused intermittent
# "could not find directory listing" failures against this mirror before.
STREAM_JOBS=$(( THREADS > 4 ? 4 : THREADS ))

BASE_URL="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp"
PANEL_URL="${BASE_URL}/release/20130502/integrated_call_samples_v3.20130502.ALL.panel"
REF_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"
REGIONS_FILE_URL="https://raw.githubusercontent.com/caterer-z-t/GRiD/main/files/734_possible_coding_vntr_regions.IBD2R_gt_0.25.uniq.txt"
REPEAT_MASK_URL="https://raw.githubusercontent.com/alexliyihao/vntrwrap/main/normalize_mosdepth/external_source/repeat_mask_list.hg38.ucsc_bed"

timestamp() { date '+%Y-%m-%d %H:%M:%S'; }
log() { echo "[$(timestamp)] $*" | tee -a "$LOG_DIR/prep.log"; }

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

log "=== Phase 0: static inputs and locus coordinates ==="

if [[ ! -f "$DATA_DIR/regions.txt" ]]; then
    log "Downloading VNTR regions file..."
    download_with_retry "$REGIONS_FILE_URL" "$DATA_DIR/regions.txt"
fi

read -r CHR START END < <(awk '$7=="LPA" {print $1, $2, $3; exit}' "$DATA_DIR/regions.txt") || true
if [[ -z "${CHR:-}" || -z "${START:-}" || -z "${END:-}" ]]; then
    echo "ERROR: Could not parse LPA coordinates from regions.txt"
    exit 1
fi
REGION="chr${CHR}:${START}-${END}"
FOCAL_BP=$(( (START + END) / 2 ))

# Persist derived coordinates so stages 2 and 3 don't have to re-derive
# (and can't silently disagree with) these values.
cat > "$DATA_DIR/pipeline_vars.sh" << VARS
export CHR="$CHR"
export START="$START"
export END="$END"
export REGION="$REGION"
export FOCAL_BP="$FOCAL_BP"
VARS

log "Locus: LPA KIV-2 ($REGION, hg38, focal_bp=$FOCAL_BP)"
log "Wrote shared coordinates to $DATA_DIR/pipeline_vars.sh"

REPEAT_MASK="$DATA_DIR/repeat_mask.bed"
if [[ ! -f "$REPEAT_MASK" ]]; then
    log "Downloading repeat mask file..."
    download_with_retry "$REPEAT_MASK_URL" "$REPEAT_MASK"
fi

log "=== Phase 1: reference genome ==="
REF_FA="$DATA_DIR/GRCh38_no_alt.fa"
if [[ ! -f "$REF_FA" ]]; then
    log "Downloading GRCh38 reference genome..."
    download_with_retry "$REF_URL" "$DATA_DIR/GRCh38_no_alt.fa.gz"
    gunzip "$DATA_DIR/GRCh38_no_alt.fa.gz"
    samtools faidx "$REF_FA"
    log "Reference genome ready: $REF_FA"
else
    log "Reference genome already present — skipping download."
fi

log "=== Phase 2: sample list ==="
PANEL="$DATA_DIR/1000G_panel.txt"
SAMPLES_FILE="$DATA_DIR/streaming_manifest.txt"
if [[ ! -f "$PANEL" ]]; then
    log "Downloading 1000G panel file..."
    download_with_retry "$PANEL_URL" "$PANEL"
fi
awk 'NR>1 {print $1, $2}' "$PANEL" > "$SAMPLES_FILE"
TOTAL=$(wc -l < "$SAMPLES_FILE")
log "Streaming manifest built: $TOTAL samples"

log "=== Phase 3: stream LPA-region CRAMs (high-coverage, 30x) ==="
FAILED_LOG="$LOG_DIR/failed_samples.txt"
: > "$FAILED_LOG"

stream_cram() {
    local sample="$1" pop="$2" ref="$3" region="$4" out_dir="$5"
    local out="${out_dir}/${sample}.cram"

    if [[ -f "$out" && -f "${out}.crai" ]]; then
        return 0
    fi

    local dir_url="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/data/${pop}/${sample}/alignment/"

    local listing="" attempt
    for attempt in 1 2 3 4 5; do
        listing=$(wget -qO- --tries=1 --timeout=30 "$dir_url" 2>/dev/null || true)
        [[ -n "$listing" ]] && break
        sleep $(( attempt * 3 + RANDOM % 3 ))
    done

    if [[ -z "$listing" ]]; then
        echo "ERROR: failed to fetch directory listing for $sample after 5 attempts: $dir_url" >&2
        echo "$sample" >> "$FAILED_LOG"
        return 1
    fi

    local full_filename
    full_filename=$(echo "$listing" | grep -oE "${sample}\.alt_bwamem_GRCh38DH\.[0-9]+\.${pop}\.low_coverage\.cram" | head -n 1)

    if [[ -z "$full_filename" ]]; then
        echo "ERROR: sample $sample has a listing but no matching high_coverage CRAM at $dir_url" >&2
        echo "$sample" >> "$FAILED_LOG"
        return 1
    fi

    local url="${dir_url}${full_filename}"

    # Each worker gets its own scratch directory so:
    #  (a) htslib's remote-index caching (which can drop a .crai named after
    #      the remote URL into the current working directory) can't collide
    #      with other concurrent workers or pollute $WORK_DIR, and
    #  (b) a stalled/hung network connection is bounded by an explicit
    #      timeout instead of silently occupying a worker slot forever —
    #      the previous version had NO timeout on this step (only on the
    #      earlier directory-listing wget), which is the likely reason most
    #      samples in this run never completed AND never got logged as
    #      failed: a stuck pipe neither succeeds nor errors.
    local scratch
    scratch=$(mktemp -d "${out_dir}/.scratch_${sample}.XXXXXX")

    (
        cd "$scratch" &&
        timeout 600 bash -c "samtools view -T '$ref' -b '$url' '$region' | samtools sort -@ 2 -o '${sample}.cram'" &&
        samtools index "${sample}.cram"
    )
    local status=$?

    if [[ $status -eq 0 && -f "$scratch/${sample}.cram" && -f "$scratch/${sample}.cram.crai" ]]; then
        mv "$scratch/${sample}.cram" "$out"
        mv "$scratch/${sample}.cram.crai" "${out}.crai"
        rm -rf "$scratch"
        return 0
    fi

    rm -rf "$scratch"
    if [[ $status -eq 124 ]]; then
        echo "ERROR: streaming timed out (600s) for $sample" >&2
    else
        echo "ERROR: streaming/sort/index failed for $sample (exit $status)" >&2
    fi
    echo "$sample" >> "$FAILED_LOG"
    return 1
}
export -f stream_cram

xargs -L 1 -P "$STREAM_JOBS" bash -c '
    ref="$1"; region="$2"; out_dir="$3"; sample="$4"; pop="$5"
    stream_cram "$sample" "$pop" "$ref" "$region" "$out_dir"
' _ "$REF_FA" "$REGION" "$CRAM_DIR" < "$SAMPLES_FILE" >> "$LOG_DIR/stream_crams.log" 2>&1 || true

N_ACTUAL=$(find "$CRAM_DIR" -maxdepth 1 -name "*.cram" | wc -l | tr -d ' ')
N_FAILED=$(wc -l < "$FAILED_LOG" | tr -d ' ')
N_UNACCOUNTED=$(( TOTAL - N_ACTUAL - N_FAILED ))
log "Streaming complete: $N_ACTUAL CRAM files in $CRAM_DIR, $N_FAILED logged failures, $N_UNACCOUNTED unaccounted for"
if [[ "$N_UNACCOUNTED" -ne 0 ]]; then
    log "WARNING: unaccounted-for samples means something neither completed nor logged a failure — check $LOG_DIR/stream_crams.log"
fi

# Clean up any stray remote-index cache files htslib may have dropped in the
# working directory during remote CRAM streaming.
find "$WORK_DIR" -maxdepth 1 \( -name "*.crai" -o -name "*.tbi" \) -delete

log "=== Phase 4: build GRiD samples file (one sample ID per line) ==="
GRID_SAMPLES_FILE="$DATA_DIR/grid_samples.txt"
: > "$GRID_SAMPLES_FILE"
for cram_path in "$CRAM_DIR"/*.cram; do
    [[ -f "$cram_path" ]] || continue
    filename=$(basename "$cram_path")
    echo "${filename%.cram}" >> "$GRID_SAMPLES_FILE"
done

N_GRID_SAMPLES=$(wc -l < "$GRID_SAMPLES_FILE" | tr -d ' ')
if [[ "$N_GRID_SAMPLES" -eq 0 ]]; then
    echo "ERROR: no CRAM files found in $CRAM_DIR — nothing to pass to GRiD."
    exit 1
fi
log "GRiD samples file ready: $GRID_SAMPLES_FILE ($N_GRID_SAMPLES samples)"

log "=== Prep complete ==="
log "Next: submit 02_run_ibs.sh, then 03_generate_config_and_run.sh"