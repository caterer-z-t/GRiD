#!/bin/bash
#SBATCH --job-name=normalize_mosdepth
#SBATCH --output=slurm/%x.out
#SBATCH --error=slurm/%x.err
#SBATCH --time=2-00:00:00
#SBATCH --partition=general
#SBATCH -n 1
#SBATCH --cpus-per-task=2
#SBATCH --mem=20G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ztcaterer@colorado.edu

# -------------------------------
# Modules
# -------------------------------
module purge

# Load Anaconda and activate environment
module load anaconda
source $(conda info --base)/etc/profile.d/conda.sh 
eval "$(conda shell.bash hook)" 
conda activate lpa_vntr

# -------------------------------
# Paths & Settings
# -------------------------------
PYTHON_SCRIPT=$LPA/playground/cram_crai/step4_normalize_mosdepth/normalize_mosdepth.py
OUT_DIR=$WORK/lpa_normalized
MOSDEPTH_DIR=$WORK/lpa_mosdepth
REPEAT_MASK_DIR=$LPA/playground/cram_crai/step4_normalize_mosdepth
N_SAMPLES=1000
OUT_FILE=$OUT_DIR/LPA_ID_scale_zdepths.txt.gz

mkdir -p $OUT_DIR

# -------------------------------
# Download repeat mask files if not present
# -------------------------------
echo "Checking for repeat mask files..."

# You'll need to download these files to longleaf. Here are some options:
# For GRCh38 (since your files are b38):
REPEAT_MASK=$REPEAT_MASK_DIR/repeat_mask_list.hg38.ucsc_bed

# Alternative options based on what you have:
# REPEAT_MASK=$REPEAT_MASK_DIR/GRCh37_repeat_mask_chr_6_cleaned.ucsc_bed  # if using GRCh37
# REPEAT_MASK=$REPEAT_MASK_DIR/repeat_mask_list.hg37_chr6.ucsc_bed         # hg37 chr6 only

if [ ! -f "$REPEAT_MASK" ]; then
    echo "ERROR: Repeat mask file not found: $REPEAT_MASK"
    echo "Please download the repeat mask files to $REPEAT_MASK_DIR"
    echo "Available options based on your mention:"
    echo "  - GRCh37_repeat_mask_chr_6_cleaned.ucsc_bed"
    echo "  - repeat_mask_list.hg37.ucsc_bed"
    echo "  - repeat_mask_list.hg37_chr6.ucsc_bed"
    echo "  - repeat_mask_list.hg38.ucsc_bed (recommended for b38 data)"
    exit 1
fi

# -------------------------------
# Count available samples
# -------------------------------
SAMPLE_COUNT=$(ls -1 $MOSDEPTH_DIR/*.regions.bed.gz | wc -l)
echo "Found $SAMPLE_COUNT .regions.bed.gz files in $MOSDEPTH_DIR"

if [ $SAMPLE_COUNT -eq 0 ]; then
    echo "ERROR: No .regions.bed.gz files found in $MOSDEPTH_DIR"
    exit 1
fi

# Use actual sample count if it's less than our max
if [ $SAMPLE_COUNT -lt $N_SAMPLES ]; then
    N_SAMPLES=$SAMPLE_COUNT
fi

echo "Processing up to $N_SAMPLES samples"

# -------------------------------
# Debugging Files
# -------------------------------
cd $OUT_DIR

# For debugging, you can run:
DEBUG_LOG=$OUT_DIR/debug_mosdepth.log
REGIONS_FILE=$LPA/playground/cram_crai/step2_count_read_files/734_possible_coding_vntr_regions.IBD2R_gt_0.25.uniq.txt

LPA_REGION=$(awk '$7=="LPA" {print $1,$2,$3,$4,$5}' "$REGIONS_FILE")
CHR=$(echo "$LPA_REGION" | awk '{print $1}')
START_HG38=$(echo "$LPA_REGION" | awk '{print $2}')
END_HG38=$(echo "$LPA_REGION" | awk '{print $3}')
START_SAMTOOLS=$(echo "$LPA_REGION" | awk '{print $4}')
END_SAMTOOLS=$(echo "$LPA_REGION" | awk '{print $5}')

echo "LPA VNTR region loaded: ${CHR}:${START_HG38}-${END_HG38}"
echo "LPA VNTR region for samtools: ${CHR}:${START_SAMTOOLS}-${END_SAMTOOLS}"

echo "Debugging mosdepth data with debug_mosdepth_data.py"
python $LPA/playground/cram_crai/step4_normalize_mosdepth/debug_mosdepth_data.py \
    $MOSDEPTH_DIR \
    $DEBUG_LOG


# -------------------------------
# Run Python Normalization
# -------------------------------

# Run the normalization script
echo "Running Python normalization..."
python $PYTHON_SCRIPT \
    $MOSDEPTH_DIR \
    $REPEAT_MASK \
    $N_SAMPLES \
    $OUT_FILE > $OUT_DIR/normalize_mosdepth_python.log 2>&1

echo "Python normalization complete. Output: $OUT_FILE"
echo "Log file: $OUT_DIR/normalize_mosdepth_python.log"

# -------------------------------
# Optional: Quick validation
# -------------------------------
if [ -f "$OUT_FILE" ]; then
    echo "Output file created successfully."
    echo "File size: $(du -h $OUT_FILE | cut -f1)"
    echo "First few lines:"
    zcat $OUT_FILE | head -5
else
    echo "ERROR: Output file not created!"
    exit 1
fi