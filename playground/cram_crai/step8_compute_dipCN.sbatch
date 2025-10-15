#!/bin/bash
#SBATCH --job-name=compute_dipCN_lpa
#SBATCH --output=slurm/%x.out
#SBATCH --error=slurm/%x.err
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --partition=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ztcaterer@colorado.edu

source "$GRID/playground/cram_crai/common.sh"

usage(){
    echo "Usage: $0 [-p PYTHON_SCRIPT] [-c COUNT_FILE] [-n NEIGHBOR_FILE] [-o OUTPUT_PREFIX] [-N N_NEIGHBORS] [-h]"
    echo ""
    echo "  -p   Python script for computing diploid CN"
    echo "  -c   Realignment count file (output from step_lpa2_realign)"
    echo "  -n   Neighbor file (output from step5_find_neighbors)"
    echo "  -o   Output prefix for diploid CN files"
    echo "  -N   Number of neighbors to use (default: 200)"
    echo "  -h   Show this help message"
    exit 1
}

# Default values
N_NEIGHBORS=200

# Parse command line arguments
while getopts "p:c:n:o:N:h" opt; do
    case $opt in
        p) COMPUTE_DIPCN_PY="$OPTARG" ;;
        c) COUNT_FILE="$OPTARG" ;;
        n) NEIGHBOR_FILE="$OPTARG" ;;
        o) OUTPUT_PREFIX="$OPTARG" ;;
        N) N_NEIGHBORS="$OPTARG" ;;
        h) usage ;;
        *) usage ;;
    esac
done

require_args COMPUTE_DIPCN_PY COUNT_FILE NEIGHBOR_FILE OUTPUT_PREFIX

export COMPUTE_DIPCN_PY COUNT_FILE NEIGHBOR_FILE OUTPUT_PREFIX N_NEIGHBORS

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
require_file "$COMPUTE_DIPCN_PY"
require_file "$COUNT_FILE"
require_file "$NEIGHBOR_FILE"

# Create output directory if needed
OUTPUT_DIR=$(dirname "$OUTPUT_PREFIX")
mkdir -p "$OUTPUT_DIR"

# Log file
LOG_FILE="${OUTPUT_DIR}/compute_dipCN.log"

# -------------------------------
# Run Analysis
# -------------------------------
echo "Starting diploid copy number computation..." > $LOG_FILE
echo "Count file: $COUNT_FILE" >> $LOG_FILE
echo "Neighbor file: $NEIGHBOR_FILE" >> $LOG_FILE
echo "Output prefix: $OUTPUT_PREFIX" >> $LOG_FILE
echo "Number of neighbors: $N_NEIGHBORS" >> $LOG_FILE
echo "" >> $LOG_FILE

start_time=$(date +%s)

# Run the Python script
python $COMPUTE_DIPCN_PY \
    --count_file "$COUNT_FILE" \
    --neighbor_file "$NEIGHBOR_FILE" \
    --output_prefix "$OUTPUT_PREFIX" \
    --n_neighbors "$N_NEIGHBORS" \
    >> $LOG_FILE 2>&1

if [ $? -eq 0 ]; then
    end_time=$(date +%s)
    elapsed=$((end_time - start_time))
    
    echo "" >> $LOG_FILE
    echo "Diploid copy number computation completed successfully!" >> $LOG_FILE
    echo "Time elapsed: ${elapsed} seconds" >> $LOG_FILE
    echo "" >> $LOG_FILE
    
    # Show summary of output files
    echo "Output files created:" >> $LOG_FILE
    for EXON in 1B_KIV3 1B_notKIV3 1B 1A; do
        OUTPUT_FILE="${OUTPUT_PREFIX}.exon${EXON}.dipCN.txt"
        if [ -f "$OUTPUT_FILE" ]; then
            NUM_LINES=$(wc -l < "$OUTPUT_FILE")
            echo "  $OUTPUT_FILE: $((NUM_LINES - 1)) samples" >> $LOG_FILE
            echo "" >> $LOG_FILE
            echo "First 10 lines of ${EXON}:" >> $LOG_FILE
            head -n 10 "$OUTPUT_FILE" >> $LOG_FILE
            echo "" >> $LOG_FILE
        fi
    done
else
    echo "ERROR: Diploid copy number computation failed!" >> $LOG_FILE
    exit 1
fi

echo "Analysis complete. Check $LOG_FILE for details."