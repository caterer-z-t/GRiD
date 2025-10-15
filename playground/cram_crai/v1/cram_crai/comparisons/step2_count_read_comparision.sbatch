#!/bin/bash
#SBATCH --job-name=count_reads_vntr_comparison
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
REGIONS_FILE=$LPA/playground/cram_crai/step2_count_read_files/734_possible_coding_vntr_regions.IBD2R_gt_0.25.uniq.txt
OUTPUT_FILE=$WORK/lpa_read_counts_comparsion/VNTR_counts.txt
REF_FASTA=$LPA/cram/hg38.fa

mkdir -p $(dirname $OUTPUT_FILE)
cd $(dirname $OUTPUT_FILE)

# -------------------------------
# Extract LPA row
# -------------------------------
LPA_LINE=$(awk '$7=="LPA"' "$REGIONS_FILE")
if [ -z "$LPA_LINE" ]; then
    echo "Error: No LPA region found!"
    exit 1
fi

CHR=$(echo "$LPA_LINE" | awk '{print $1}')
HG38_START=$(echo "$LPA_LINE" | awk '{print $2}')
HG38_END=$(echo "$LPA_LINE" | awk '{print $3}')
SAM_START=$(echo "$LPA_LINE" | awk '{print $4}')
SAM_END=$(echo "$LPA_LINE" | awk '{print $5}')

echo "HG38 region: ${CHR}:${HG38_START}-${HG38_END}"
echo "Samtools region: ${CHR}:${SAM_START}-${SAM_END}"

# -------------------------------
# Start Counting
# -------------------------------
start_time=$(date +%s)

CRAM_FILES=($CRAM_DIR/*.cram)
echo "Detected ${#CRAM_FILES[@]} CRAM files."

# Header
echo -e "Sample\tHG38_Count\tSamtools_Count\tDifference\tPercent_Diff" > "$OUTPUT_FILE"

# Loop over CRAMs
for CRAM in "${CRAM_FILES[@]}"; do
    BASENAME=$(basename "$CRAM")
    echo "Processing $BASENAME..."

    # Ensure index exists
    if [ ! -f "${CRAM}.crai" ]; then
        echo "Index not found for $BASENAME. Creating..."
        samtools index -c "$CRAM"
    fi

    # Count properly paired reads in HG38 region
    HG38_COUNT=$(samtools view -T "$REF_FASTA" "$CRAM" "chr${CHR}:${HG38_START}-${HG38_END}" \
        | awk '$2==83 || $2==99 || $2==147 || $2==163 || $2==81 || $2==97 || $2==145 || $2==161 {ctr++} END {print ctr+0}')

    # Count properly paired reads in Samtools region
    SAM_COUNT=$(samtools view -T "$REF_FASTA" "$CRAM" "chr${CHR}:${SAM_START}-${SAM_END}" \
        | awk '$2==83 || $2==99 || $2==147 || $2==163 || $2==81 || $2==97 || $2==145 || $2==161 {ctr++} END {print ctr+0}')

    # Difference and percent diff
    DIFF=$((HG38_COUNT - SAM_COUNT))
    if [ "$HG38_COUNT" -ne 0 ]; then
        PCT_DIFF=$(awk -v d="$DIFF" -v h="$HG38_COUNT" 'BEGIN {printf "%.2f", (d/h*100)}')
    else
        PCT_DIFF="NA"
    fi

    echo -e "${BASENAME}\t${HG38_COUNT}\t${SAM_COUNT}\t${DIFF}\t${PCT_DIFF}" >> "$OUTPUT_FILE"
done

echo "VNTR read counting completed in $(($(date +%s)-$start_time)) seconds."

# -------------------------------
# Mosdepth Comparison
# -------------------------------
plot_py="$LPA/playground/cram_crai/step2_count_read_files/plot_dist.py"
python $plot_py VNTR_counts.txt -o VNTR
