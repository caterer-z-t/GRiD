#!/bin/bash

module load anaconda mosdepth samtools 

POP_FILTER=""
N_SAMPLES=0   # 0 = all
JOBS_OVERRIDE=0
 
usage() {
    grep '^#' "$0" | sed -e 's/^#//' -e '1d'
    exit 0
}
 
while [[ $# -gt 0 ]]; do
    case $1 in
        --pop)  POP_FILTER="$2"; shift 2 ;;
        --n)    N_SAMPLES="$2";  shift 2 ;;
        --jobs) JOBS_OVERRIDE="$2"; shift 2 ;;
        -h|--help) usage ;;
        *) echo "Unknown argument: $1"; usage ;;
    esac
done
 
# Dependency check
MISSING=()
for cmd in samtools mosdepth conda awk wget; do
    command -v "$cmd" &>/dev/null || MISSING+=("$cmd")
done
if [[ ${#MISSING[@]} -gt 0 ]]; then
    echo "ERROR: missing required tools: ${MISSING[*]}"
    echo "Install with: conda install -c bioconda ${MISSING[*]}"
    exit 1
fi
 
# Paths and parameters
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
IBS_SCRIPT="$SCRIPT_DIR/IBS_example.sh"
 
WORK_DIR="$(pwd)"
CRAM_DIR="$WORK_DIR/crams"
LOG_DIR="$WORK_DIR/logs"
DATA_DIR="$WORK_DIR/data"
OUTPUT_DIR="$WORK_DIR/output"
MOSDEPTH_WORK="$WORK_DIR/mosdepth_work"
 
mkdir -p "$CRAM_DIR" "$LOG_DIR" "$DATA_DIR" "$OUTPUT_DIR" "$MOSDEPTH_WORK"
 
THREADS="${SLURM_CPUS_PER_TASK:-$(nproc)}"
if [[ "$JOBS_OVERRIDE" -gt 0 ]]; then
    STREAM_JOBS="$JOBS_OVERRIDE"
else
    STREAM_JOBS=$(( THREADS > 2 ? 2 : THREADS ))
fi
 
BASE_URL="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp"
PANEL_URL="${BASE_URL}/release/20130502/integrated_call_samples_v3.20130502.ALL.panel"
REF_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"
REGIONS_FILE_URL="https://raw.githubusercontent.com/caterer-z-t/GRiD/main/files/734_possible_coding_vntr_regions.IBD2R_gt_0.25.uniq.txt"
REPEAT_MASK_URL="https://raw.githubusercontent.com/alexliyihao/vntrwrap/main/normalize_mosdepth/external_source/repeat_mask_list.hg38.ucsc_bed"
 
timestamp() { date '+%Y-%m-%d %H:%M:%S'; }
log() { echo "[$(timestamp)] $*" | tee -a "$LOG_DIR/pipeline.log"; }
 
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
 
# Fetch static inputs
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

REPEAT_MASK="$DATA_DIR/repeat_mask.bed"

if [[ ! -f "$REPEAT_MASK" ]]; then
    log "Downloading repeat mask file..."
    download_with_retry "$REPEAT_MASK_URL" "$REPEAT_MASK"
fi
 
log "  GRiD — 1000 Genomes Example"
log "  Locus:   LPA KIV-2  ($REGION, hg38)"
log "  Threads: $THREADS  |  Parallel streams: $STREAM_JOBS"
log "  Output:  $OUTPUT_DIR"
 
# Activate environment
if ! conda env list | grep -qE '^\s*grid\s'; then
    echo "ERROR: conda environment 'grid' not found."
    echo "Create it first, e.g.: conda env create -f environment.yml"
    exit 1
fi
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate grid
 
# Phase 1 — Reference genome
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
 
# Phase 2 — Sample list
PANEL="$DATA_DIR/1000G_panel.txt"
SAMPLES_FILE="$DATA_DIR/streaming_manifest.txt"
 
if [[ ! -f "$PANEL" ]]; then
    log "Downloading 1000G panel file..."
    download_with_retry "$PANEL_URL" "$PANEL"
fi
 
if [[ -n "$POP_FILTER" ]]; then
    log "Filtering to superpopulation: $POP_FILTER"
    awk -v pop="$POP_FILTER" 'NR>1 && $3==pop {print $1, $2}' "$PANEL" > "$SAMPLES_FILE"
else
    awk 'NR>1 {print $1, $2}' "$PANEL" > "$SAMPLES_FILE"
fi
 
TOTAL=$(wc -l < "$SAMPLES_FILE")
if [[ "$N_SAMPLES" -gt 0 && "$N_SAMPLES" -lt "$TOTAL" ]]; then
    head -n "$N_SAMPLES" "$SAMPLES_FILE" > "${SAMPLES_FILE}.tmp" && mv "${SAMPLES_FILE}.tmp" "$SAMPLES_FILE"
    log "Using first $N_SAMPLES of $TOTAL available samples."
else
    log "Using all $TOTAL samples."
fi
 
# Phase 3 — Stream LPA-region CRAMs (parallel, remote)
log "Streaming LPA region CRAMs from 1000G EBI..."
 
FAILED_LOG="$LOG_DIR/failed_samples.txt"
: > "$FAILED_LOG"
 
stream_cram() {
    local sample="$1"
    local pop="$2"
    local ref="$3"
    local region="$4"
    local out_dir="$5"
 
    local out="${out_dir}/${sample}.cram"
 
    if [[ -f "$out" && -f "${out}.crai" ]]; then
        return 0
    fi
 
    local dir_url="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/data/${pop}/${sample}/alignment/"
 
    local listing=""
    local attempt
    for attempt in 1 2 3 4 5; do
        listing=$(wget -qO- --tries=1 --timeout=30 "$dir_url" 2>/dev/null || true)
        if [[ -n "$listing" ]]; then
            break
        fi
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
        echo "ERROR: sample $sample has a listing but no matching low_coverage CRAM at $dir_url" >&2
        echo "$sample" >> "$FAILED_LOG"
        return 1
    fi
 
    local url="${dir_url}${full_filename}"
 
    if ! (samtools view -T "$ref" -b "$url" "$region" | \
          samtools sort -@ 2 -o "$out" && \
          samtools index "$out"); then
        echo "ERROR: streaming/sort/index failed for $sample — cleaning up partial output" >&2
        rm -f "$out" "${out}.crai"
        echo "$sample" >> "$FAILED_LOG"
        return 1
    fi
}
 
export -f stream_cram
xargs -L 1 -P "$STREAM_JOBS" bash -c '
    ref="$1"
    region="$2"
    out_dir="$3"
    sample="$4"
    pop="$5"
    stream_cram "$sample" "$pop" "$ref" "$region" "$out_dir"
' _ "$REF_FA" "$REGION" "$CRAM_DIR" < "$SAMPLES_FILE" >> "$LOG_DIR/stream_crams.log" 2>&1 || true

find "$WORK_DIR" -maxdepth 1 -name "*.low_coverage.cram.crai" -delete

# Compile GRiD Input manifest file directly from downloaded assets
GRID_SAMPLES_FILE="$DATA_DIR/grid_samples.txt"
: > "$GRID_SAMPLES_FILE"

for cram_path in "$CRAM_DIR"/*.cram; do
    if [[ -f "$cram_path" ]]; then
        filename=$(basename "$cram_path")
        sample_id="${filename%.cram}"
        echo "$sample_id" >> "$GRID_SAMPLES_FILE"
    fi
done

N_GRID_SAMPLES=$(wc -l < "$GRID_SAMPLES_FILE" | tr -d ' ')
if [[ "$N_GRID_SAMPLES" -eq 0 ]]; then
    echo "ERROR: no files found in $CRAM_DIR — nothing to pass to GRiD."
    exit 1
fi

# Phase 4 — Run IBS_example.sh if needed
IBS_OUTPUT="$DATA_DIR/ibs_neighbors_chr6.tsv.gz"
HAPLOID_RUN="False"
 
if [[ -f "$IBS_OUTPUT" ]]; then
    log "IBS neighbors file already present — skipping IBS_example.sh."
    HAPLOID_RUN="True"
else
    if [[ ! -f "$IBS_SCRIPT" ]]; then
        log "WARNING: IBS neighbors file not found and $IBS_SCRIPT is missing."
        log "         Haploid CN estimation will be skipped for this run."
    else
        log "IBS neighbors file not found — running IBS_example.sh to generate it..."
        if WORK_DIR="$WORK_DIR" bash "$IBS_SCRIPT" --focal-bp "$FOCAL_BP"; then
            if [[ -f "$IBS_OUTPUT" ]]; then
                log "IBS neighbors file generated: $IBS_OUTPUT"
                HAPLOID_RUN="True"
            else
                log "WARNING: IBS_example.sh completed but $IBS_OUTPUT still missing."
                log "         Haploid CN estimation will be skipped for this run."
            fi
        else
            log "WARNING: IBS_example.sh failed. Haploid CN estimation will be skipped for this run."
            log "         See $LOG_DIR/IBS.log for details."
        fi
    fi
fi

# Phase 5 — Generate config
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
  run: $HAPLOID_RUN
  output_file_prefix: "haploid_genotypes"
  ibs_output: "$DATA_DIR/ibs_neighbors_chr6.tsv.gz"
  min_neighbors: 1
  max_neighbors: 10
  n_iters: 100
YAML

log "Config written to: $CONFIG"

# Phase 5 — Run GRiD
log "Running GRiD pipeline..."
grid wgs $CONFIG