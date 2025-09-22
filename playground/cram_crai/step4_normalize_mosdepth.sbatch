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

source "$(dirname "$0")/common.sh"

usage() {
    echo "Usage: $0 [-p PYTHON_SCRIPT] [-f OUTPUT_FILE] [-d OUTPUT_DIR] [-m MOSDEPTH_PREFIX] [-r REPEAT_MASK] [-e EXAMPLE_REGIONS] [-n N_SAMPLES] [-h]"
    echo ""
    echo "  -p   Python script to run"
    echo "  -f   Output file path"
    echo "  -d   Output directory"
    echo "  -m   Mosdepth prefix (for batch files: prefix_batch_N.txt.gz)"
    echo "  -r   Repeat mask BED file"
    echo "  -e   Example regions file (.regions.bed.gz)"
    echo "  -n   Number of samples"
    echo "  -h   Show this help message"
    exit 1
}

# Parse command line arguments
while getopts "p:f:d:m:r:e:n:h" opt; do
    case $opt in
        p) PYTHON_SCRIPT="$OPTARG" ;;
        f) OUT_FILE="$OPTARG" ;;
        d) OUT_DIR="$OPTARG" ;;
        m) MOSDEPTH_PREFIX="$OPTARG" ;;
        r) REPEAT_MASK="$OPTARG" ;;
        e) EXAMPLE_REGIONS="$OPTARG" ;;
        n) N_SAMPLES="$OPTARG" ;;
        h) usage ;;
        *) usage ;;
    esac
done

require_args PYTHON_SCRIPT OUT_FILE OUT_DIR MOSDEPTH_PREFIX REPEAT_MASK EXAMPLE_REGIONS N_SAMPLES

export PYTHON_SCRIPT OUT_FILE OUT_DIR MOSDEPTH_PREFIX REPEAT_MASK EXAMPLE_REGIONS N_SAMPLES

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
require_file "$PYTHON_SCRIPT"
require_file "$REPEAT_MASK"
require_file "$EXAMPLE_REGIONS"

require_dir "$OUT_DIR"

# Check that at least one batch file exists
if [ ! -f "${MOSDEPTH_PREFIX}_batch_1.txt.gz" ]; then
    echo "Error: No batch files found with prefix '$MOSDEPTH_PREFIX'"
    echo "Expected files like: ${MOSDEPTH_PREFIX}_batch_1.txt.gz"
    exit 1
fi

mkdir -p $OUT_DIR

# -------------------------------
# Run Python Normalization
# -------------------------------

# Run the normalization script with correct arguments
echo "Running Python normalization..."
python $PYTHON_SCRIPT \
    $MOSDEPTH_PREFIX \
    $REPEAT_MASK \
    $EXAMPLE_REGIONS \
    $N_SAMPLES \
    $OUT_FILE > $OUT_DIR/normalize_mosdepth_python.log 2>&1

if [ $? -eq 0 ]; then
    echo "Python normalization complete. Output: $OUT_FILE"
    echo "Log file: $OUT_DIR/normalize_mosdepth_python.log"
else
    echo "Error: Python normalization failed. Check log file: $OUT_DIR/normalize_mosdepth_python.log"
    exit 1
fi