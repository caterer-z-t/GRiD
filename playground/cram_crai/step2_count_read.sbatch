#!/bin/bash
#SBATCH --job-name=count_reads_vntr
#SBATCH --output=slurm/%x.out
#SBATCH --error=slurm/%x.err
#SBATCH --time=2-00:00:00
#SBATCH --partition=general
#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ztcaterer@colorado.edu

pwd

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
# Modules -- maybe change to singularity or docker env?
# -------------------------------
module purge
module load samtools/1.22

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
# Extract only LPA row
# -------------------------------
CHR=$(normalize_chr "$CHR")

echo "Counting reads for LPA VNTR: ${CHR}:${START}-${END}"

# -------------------------------
# Start Counting
# -------------------------------
start_time=$(date +%s)

# Get CRAM files from directory
CRAM_FILES=($(find_cram_files "$CRAM_DIR"))

echo "Detected ${#CRAM_FILES[@]} CRAM files."

# Write header: Sample + LPA region
echo -e "Sample\t${CHR}:${START}-${END}" > "$OUTPUT_FILE"

# Loop over each CRAM
for CRAM in "${CRAM_FILES[@]}"; do
    BASENAME=$(basename "$CRAM")

    echo "Processing $BASENAME..."

    # Create CRAM index if missing
    if [ ! -f "${CRAM}.crai" ]; then
        echo "Index not found for $BASENAME. Creating..."
        samtools index -c "$CRAM"
    fi

    # Count properly paired reads in LPA VNTR
    COUNT=$(samtools view -T "$REF_FASTA" "$CRAM" "${CHR}:${START}-${END}" \
            | awk '$2==83 || $2==99 || $2==147 || $2==163 || $2==81 || $2==97 || $2==145 || $2==161 {ctr++} END {print ctr+0}')

    echo -e "${BASENAME}\t${COUNT}" >> "$OUTPUT_FILE"
done

echo "VNTR read counting completed in $(($(date +%s)-$start_time)) seconds."