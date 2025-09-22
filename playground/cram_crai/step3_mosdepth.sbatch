#!/bin/bash
#SBATCH --job-name=mosdepth_lpa
#SBATCH --output=slurm/%x.out
#SBATCH --error=slurm/%x.err
#SBATCH --time=2-00:00:00
#SBATCH --partition=general
#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ztcaterer@colorado.edu

source "$LPA/playground/cram_crai/common.sh"

# Parse command line arguments
while getopts "c:o:r:C:s:e:w:h" opt; do
    case $opt in
        c) CRAM_DIR="$OPTARG" ;;
        o) OUTPUT_FILE="$OPTARG" ;;
        r) REF_FASTA="$OPTARG" ;;
        C) CHR="$OPTARG" ;;
        s) START="$OPTARG" ;;
        e) END="$OPTARG" ;;
        w) WORK_DIR="$OPTARG" ;;
        h) usage_common ;;
        *) usage_common ;;
    esac
done

# Run checks via common.sh helpers
require_args CRAM_DIR OUTPUT_FILE REF_FASTA CHR START END
require_dir "$CRAM_DIR"
require_file "$REF_FASTA"

export CRAM_DIR OUTPUT_FILE REF_FASTA CHR START END WORK_DIR

# -------------------------------
# Modules & Paths
# -------------------------------
module purge
module load mosdepth/0.3.11

# -------------------------------
# Paths & Settings
# -------------------------------

# Normalize WORK_DIR relative to OUTPUT_FILE
WORK_DIR=$(sync_workdir_with_output "$OUTPUT_FILE" "$WORK_DIR")

# Create working directory if it doesn't exist
mkdir -p "$WORK_DIR"

# Move to working directory
cd "$WORK_DIR"

# -------------------------------
# Extract LPA VNTR coordinates
# -------------------------------

CHR=$(normalize_chr "$CHR")

echo "Counting reads for LPA VNTR: ${CHR}:${START}-${END}"

# -------------------------------
# Start Mosdepth
# -------------------------------
start_time=$(date +%s)

CRAM_FILES=($(find_cram_files "$CRAM_DIR"))

echo "Detected ${#CRAM_FILES[@]} CRAM files."

# Write header: Sample + LPA region
echo -e "Sample\t${CHR}:${START}-${END}" > "$OUTPUT_FILE"

for CRAM in "${CRAM_FILES[@]}"; do
    BASENAME=$(basename $CRAM .cram)
    echo "Processing $BASENAME..."

    # Ensure index exists
    if [ ! -f "${CRAM}.crai" ]; then
        echo "Index not found for $BASENAME. Creating..."
        samtools index -c "$CRAM"
    fi

    OUT_PREFIX=$WORK_DIR/${BASENAME}_LPA

    mosdepth -n --fast-mode \
        --by 1000 \
        -f $REF_FASTA \
        $OUT_PREFIX \
        $CRAM

    # Extract coverage from mosdepth output for region of interest
    REGION_COV=$(zgrep -P "^${CHR}\t${START}\t${END}" ${OUT_PREFIX}.regions.bed.gz | \
                 awk '{printf("%d", int(100*$4+0.5))}')

    echo -e "${BASENAME}\t${REGION_COV}" >> "$OUTPUT_FILE"
    
done

echo "Mosdepth VNTR coverage completed in $(($(date +%s)-$start_time)) seconds."