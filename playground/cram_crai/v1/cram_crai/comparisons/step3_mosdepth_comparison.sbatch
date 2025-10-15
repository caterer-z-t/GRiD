#!/bin/bash
#SBATCH --job-name=mosdepth_lpa_compare
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
module load samtools

module load anaconda
source $(conda info --base)/etc/profile.d/conda.sh 
eval "$(conda shell.bash hook)"
conda activate lpa_vntr

CRAM_DIR=$SOL/cram
REF_FASTA=$LPA/cram/hg38.fa
OUTPUT_DIR=$WORK/lpa_mosdepth_comparison
REGIONS_FILE="$LPA/playground/cram_crai/step2_count_read_files/734_possible_coding_vntr_regions.IBD2R_gt_0.25.uniq.txt"
COMPARE_SCRIPT="$LPA/playground/cram_crai/step3_mosdepth_comparison/compare_mosdepth_cov.py"
LOG_FILE=$OUTPUT_DIR/mosdepth_lpa_compare.log

mkdir -p $OUTPUT_DIR
cd $OUTPUT_DIR

# -------------------------------
# Extract LPA row
# -------------------------------
echo "Extracting LPA region from $REGIONS_FILE ..."
LPA_LINE=$(awk '$7=="LPA"' "$REGIONS_FILE")

if [ -z "$LPA_LINE" ]; then
    echo "Error: No LPA region found!"
    exit 1
fi

# Parse coordinates
read CHR START1 END1 START2 END2 REST <<< "$LPA_LINE"

echo "HG38 region: ${CHR}:${START1}-${END1}"
echo "Samtools region: ${CHR}:${START2}-${END2}"

# -------------------------------
# Initialize results file
# -------------------------------
OUT_TABLE=$OUTPUT_DIR/LPA_regions_comparison.txt
echo -e "SAMPLE\tREGION1_MEAN\tREGION2_MEAN\tDIFFERENCE\tPERCENT_DIFF" > $OUT_TABLE

# -------------------------------
# Process each CRAM
# -------------------------------
for CRAM in $CRAM_DIR/*.cram; do
    BASENAME=$(basename $CRAM .cram)
    echo "Processing $BASENAME..."
    
    # Ensure index exists
    if [ ! -f "${CRAM}.crai" ]; then
        echo "Index not found for $BASENAME. Creating..."
        samtools index -c "$CRAM"
    fi
    
    OUT_PREFIX=$OUTPUT_DIR/${BASENAME}_LPA
    
    # Run mosdepth
    mosdepth -n --fast-mode \
        --by 1000 \
        -f $REF_FASTA \
        $OUT_PREFIX \
        $CRAM

    # Compare regions with Python
    python $COMPARE_SCRIPT \
        $OUT_PREFIX.regions.bed.gz \
        $CHR $START1 $END1 \
        $START2 $END2 \
        $OUT_TABLE \
        $LOG_FILE
done
