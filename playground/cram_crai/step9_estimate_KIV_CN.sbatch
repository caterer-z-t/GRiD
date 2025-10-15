#!/bin/bash
#SBATCH --job-name=estimate_kiv_cn
#SBATCH --output=slurm/%x.out
#SBATCH --error=slurm/%x.err
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --partition=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ztcaterer@colorado.edu

source "$GRID/playground/cram_crai/common.sh"

usage(){
    echo "Usage: $0 [-p PYTHON_SCRIPT] [-a EXON1A_FILE] [-b EXON1B_FILE] [-o OUTPUT_FILE] [-f FORMAT] [-h]"
    echo ""
    echo "  -p   Python script for estimating KIV CN"
    echo "  -a   Exon1A diploid CN file (*.exon1A.dipCN.txt)"
    echo "  -b   Exon1B diploid CN file (*.exon1B.dipCN.txt)"
    echo "  -o   Output file path"
    echo "  -f   Output format: csv, tsv, or txt (default: tsv)"
    echo "  -h   Show this help message"
    exit 1
}

# Default values
OUTPUT_FORMAT="tsv"

# Parse command line arguments
while getopts "p:a:b:o:f:h" opt; do
    case $opt in
        p) ESTIMATE_KIV_PY="$OPTARG" ;;
        a) EXON1A_FILE="$OPTARG" ;;
        b) EXON1B_FILE="$OPTARG" ;;
        o) OUTPUT_FILE="$OPTARG" ;;
        f) OUTPUT_FORMAT="$OPTARG" ;;
        h) usage ;;
        *) usage ;;
    esac
done

require_args ESTIMATE_KIV_PY EXON1A_FILE EXON1B_FILE OUTPUT_FILE

export ESTIMATE_KIV_PY EXON1A_FILE EXON1B_FILE OUTPUT_FILE OUTPUT_FORMAT

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
require_file "$ESTIMATE_KIV_PY"
require_file "$EXON1A_FILE"
require_file "$EXON1B_FILE"

# Create output directory if needed
OUTPUT_DIR=$(dirname "$OUTPUT_FILE")
mkdir -p "$OUTPUT_DIR"

# Log file
LOG_FILE="${OUTPUT_DIR}/estimate_KIV_CN.log"

# -------------------------------
# Run Analysis
# -------------------------------
echo "Starting KIV2 copy number estimation..." > $LOG_FILE
echo "Exon1A file: $EXON1A_FILE" >> $LOG_FILE
echo "Exon1B file: $EXON1B_FILE" >> $LOG_FILE
echo "Output file: $OUTPUT_FILE" >> $LOG_FILE
echo "Output format: $OUTPUT_FORMAT" >> $LOG_FILE
echo "" >> $LOG_FILE

start_time=$(date +%s)

# Run the Python script
python $ESTIMATE_KIV_PY \
    --exon1a "$EXON1A_FILE" \
    --exon1b "$EXON1B_FILE" \
    --output "$OUTPUT_FILE" \
    --format "$OUTPUT_FORMAT" \
    >> $LOG_FILE 2>&1

if [ $? -eq 0 ]; then
    end_time=$(date +%s)
    elapsed=$((end_time - start_time))
    
    echo "" >> $LOG_FILE
    echo "KIV2 copy number estimation completed successfully!" >> $LOG_FILE
    echo "Time elapsed: ${elapsed} seconds" >> $LOG_FILE
    echo "" >> $LOG_FILE
    
    # Show summary statistics
    NUM_SAMPLES=$(tail -n +2 "$OUTPUT_FILE" | wc -l)
    echo "Total samples processed: $NUM_SAMPLES" >> $LOG_FILE
    echo "" >> $LOG_FILE
    
    # Show first few lines of output
    echo "First 20 lines of output file:" >> $LOG_FILE
    head -n 20 "$OUTPUT_FILE" >> $LOG_FILE
    
    echo "" >> $LOG_FILE
    echo "Output file: $OUTPUT_FILE" >> $LOG_FILE
else
    echo "ERROR: KIV2 copy number estimation failed!" >> $LOG_FILE
    exit 1
fi

echo "Analysis complete. Check $LOG_FILE for details."