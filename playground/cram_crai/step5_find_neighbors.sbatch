#!/bin/bash
#SBATCH --job-name=find_neighbors_lpa
#SBATCH --output=slurm/%x.out
#SBATCH --error=slurm/%x.err
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --partition=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ztcaterer@colorado.edu

source "$LPA/playground/cram_crai/common.sh"

usage(){
    echo "Usage: $0 [-p PYTHON_SCRIPT] [-h]"
    echo ""
    echo "  -p   Python script to run"
    echo "  -o   Output directory"
    echo "  -i   Input file with normalized z-depths"
    echo "  -f   Output file prefix"
    echo "  -z   Max z-score distance to consider a neighbor"
    echo "  -n   Number of nearest neighbors to find"
    echo "  -h   Show this help message"
    exit 1
}

# Parse command line arguments
while getopts "p:o:i:f:z:n:h" opt; do
    case $opt in
        p) FIND_NEIGHBORS_PY="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        i) INPUT_FILE="$OPTARG" ;;
        f) OUTPUT_PREFIX="$OPTARG" ;;
        z) MAX_Z_SCORE="$OPTARG" ;;
        n) N_NEIGHBORS="$OPTARG" ;;
        h) usage ;;
        *) usage ;;
    esac
done

require_args FIND_NEIGHBORS_PY OUTPUT_DIR INPUT_FILE OUTPUT_PREFIX MAX_Z_SCORE N_NEIGHBORS

export FIND_NEIGHBORS_PY OUTPUT_DIR INPUT_FILE OUTPUT_PREFIX MAX_Z_SCORE N_NEIGHBORS


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
require_file "$FIND_NEIGHBORS_PY"
require_file "$INPUT_FILE"

mkdir -p $OUTPUT_DIR

cd $OUTPUT_DIR

# -------------------------------
# Run Analysis
# -------------------------------

# Run the Python script
echo "Starting find_neighbors analysis..." > $OUTPUT_DIR/find_neighbors.log 2>&1
python $FIND_NEIGHBORS_PY \
    --input_file $INPUT_FILE \
    --output_prefix $OUTPUT_PREFIX \
    --zmax $MAX_Z_SCORE \
    --n_neighbors $N_NEIGHBORS >> $OUTPUT_DIR/find_neighbors.log 2>&1

echo "Find neighbors analysis completed. Results in $OUTPUT_DIR" >> $OUTPUT_DIR/find_neighbors.log 2>&1
# print the first few lines of the output files

echo "First few lines of the output file:" >> $OUTPUT_DIR/find_neighbors.log 2>&1

zcat "${OUTPUT_PREFIX}.zMax${MAX_Z_SCORE}.txt.gz" | head -n 20 >> $OUTPUT_DIR/find_neighbors.log 2>&1