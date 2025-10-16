#!/bin/bash
#SBATCH --job-name=plot_lpa_assoc
#SBATCH --output=slurm/%x.out
#SBATCH --error=slurm/%x.err
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --partition=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ztcaterer@colorado.edu

source "$GRID/playground/cram_crai/common.sh"

usage(){
    echo "Usage: $0 [-p PYTHON_SCRIPT] [-k KIV2_FILE] [-l LPA_FILE] [-o OUTPUT_PREFIX] [-h]"
    echo ""
    echo "  -p   Python plotting script"
    echo "  -k   KIV2 copy number estimates file (output from step9)"
    echo "  -l   Lp(a) measurements file (diploid_calls.txt)"
    echo "  -o   Output prefix for plots and merged data"
    echo "  -h   Show this help message"
    exit 1
}

# Parse command line arguments
while getopts "p:k:l:o:h" opt; do
    case $opt in
        p) PLOT_PY="$OPTARG" ;;
        k) KIV2_FILE="$OPTARG" ;;
        l) LPA_FILE="$OPTARG" ;;
        o) OUTPUT_PREFIX="$OPTARG" ;;
        h) usage ;;
        *) usage ;;
    esac
done

require_args PLOT_PY KIV2_FILE LPA_FILE OUTPUT_PREFIX

export PLOT_PY KIV2_FILE LPA_FILE OUTPUT_PREFIX

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
require_file "$PLOT_PY"
require_file "$KIV2_FILE"
require_file "$LPA_FILE"

# Create output directory if needed
OUTPUT_DIR=$(dirname "$OUTPUT_PREFIX")
mkdir -p "$OUTPUT_DIR"

# Log file
LOG_FILE="${OUTPUT_DIR}/plot_association.log"

# -------------------------------
# Run Analysis
# -------------------------------
echo "Starting Lp(a) association plotting..." > $LOG_FILE
echo "KIV2 estimates: $KIV2_FILE" >> $LOG_FILE
echo "Lp(a) data: $LPA_FILE" >> $LOG_FILE
echo "Output prefix: $OUTPUT_PREFIX" >> $LOG_FILE
echo "" >> $LOG_FILE

start_time=$(date +%s)

# Run the Python script
python $PLOT_PY \
    --kiv2_file "$KIV2_FILE" \
    --lpa_file "$LPA_FILE" \
    --output_prefix "$OUTPUT_PREFIX" \
    >> $LOG_FILE 2>&1

if [ $? -eq 0 ]; then
    end_time=$(date +%s)
    elapsed=$((end_time - start_time))
    
    echo "" >> $LOG_FILE
    echo "Plotting completed successfully!" >> $LOG_FILE
    echo "Time elapsed: ${elapsed} seconds" >> $LOG_FILE
    echo "" >> $LOG_FILE
    
    # List generated files
    echo "Generated files:" >> $LOG_FILE
    ls -lh "${OUTPUT_PREFIX}"* 2>/dev/null | awk '{print "  " $9 " (" $5 ")"}' >> $LOG_FILE
    
else
    echo "ERROR: Plotting failed!" >> $LOG_FILE
    exit 1
fi

echo "Analysis complete. Check $LOG_FILE for details."
echo "Plots saved to: ${OUTPUT_PREFIX}_*.png"