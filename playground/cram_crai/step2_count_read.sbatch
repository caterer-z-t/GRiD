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

# -------------------------------
# Modules
# -------------------------------
module purge
module load samtools/1.22

# -------------------------------
# Paths & Settings
# -------------------------------
CRAM_DIR=$SOL/cram
IBD2R_PATH=$LPA/playground/cram_crai/step2_count_read_files/734_possible_coding_vntr_regions.IBD2R_gt_0.25.uniq.txt
OUTPUT_FILE=$WORK/lpa_read_counts/VNTR_counts.txt
REF_FASTA=$LPA/cram/hg38.fa

mkdir -p $(dirname $OUTPUT_FILE)

cd $(dirname $OUTPUT_FILE)

# -------------------------------
# Extract only LPA row
# -------------------------------
LPA_REGION=$(awk '$7=="LPA" {print $1,$4,$5}' "$IBD2R_PATH")
CHR=$(echo "$LPA_REGION" | awk '{print $1}')
START=$(echo "$LPA_REGION" | awk '{print $2}')
END=$(echo "$LPA_REGION" | awk '{print $3}')

echo "Counting reads for LPA VNTR: ${CHR}:${START}-${END}"

# -------------------------------
# Start Counting
# -------------------------------
start_time=$(date +%s)

# Find all CRAM files
CRAM_FILES=($CRAM_DIR/*.cram)
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
    COUNT=$(samtools view -T "$REF_FASTA" "$CRAM" "chr${CHR}:${START}-${END}" \
            | awk '$2==83 || $2==99 || $2==147 || $2==163 || $2==81 || $2==97 || $2==145 || $2==161 {ctr++} END {print ctr+0}')

    echo -e "${BASENAME}\t${COUNT}" >> "$OUTPUT_FILE"
done

echo "VNTR read counting completed in $(($(date +%s)-$start_time)) seconds."