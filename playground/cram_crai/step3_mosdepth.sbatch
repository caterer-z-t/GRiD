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

# -------------------------------
# Modules & Paths
# -------------------------------
module purge
module load mosdepth/0.3.11

CRAM_DIR=$SOL/cram
REF_FASTA=$LPA/cram/hg38.fa
OUTPUT_DIR=$WORK/lpa_mosdepth
REGIONS_FILE="$LPA/playground/cram_crai/step2_count_read_files/734_possible_coding_vntr_regions.IBD2R_gt_0.25.uniq.txt"

mkdir -p $OUTPUT_DIR
cd $OUTPUT_DIR

# -------------------------------
# Extract LPA VNTR coordinates
# -------------------------------
LPA_REGION=$(awk '$7=="LPA" {print $1,$2,$3}' "$REGIONS_FILE")
CHR=$(echo "$LPA_REGION" | awk '{print $1}')
START=$(echo "$LPA_REGION" | awk '{print $2}')
END=$(echo "$LPA_REGION" | awk '{print $3}')

echo "LPA VNTR region loaded: ${CHR}:${START}-${END}"

# -------------------------------
# Start Mosdepth
# -------------------------------
start_time=$(date +%s)

for CRAM in $CRAM_DIR/*.cram; do
    BASENAME=$(basename $CRAM .cram)
    echo "Processing $BASENAME..."

    # Ensure index exists
    if [ ! -f "${CRAM}.crai" ]; then
        echo "Index not found for $BASENAME. Creating..."
        samtools index -c "$CRAM"
    fi

    OUT_PREFIX=$OUTPUT_DIR/${BASENAME}_LPA

    mosdepth -n --fast-mode \
        --by 1000 \
        -f $REF_FASTA \
        $OUT_PREFIX \
        $CRAM

    # Extract only chr6 coverage for LPA VNTR
    zcat ${OUT_PREFIX}.regions.bed.gz \
        | awk -v start=$START -v end=$END -v chr=$CHR '$1==chr && $2>=start && $3<=end { sum+=$4 } END { print BASENAME "\t" sum }' \
        BASENAME=$BASENAME \
        >> $OUTPUT_DIR/LPA_mosdepth_summary.txt
done

echo "Mosdepth VNTR coverage completed in $(($(date +%s)-$start_time)) seconds."