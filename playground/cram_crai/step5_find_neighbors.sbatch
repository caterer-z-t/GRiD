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

# Load required modules
module purge
module load anaconda
source $(conda info --base)/etc/profile.d/conda.sh 
eval "$(conda shell.bash hook)" 

conda activate lpa_vntr

# Set variables - MODIFY THESE PATHS FOR YOUR SETUP
BATCH_NUM=0
TOTAL_BATCHES=1
MAX_Z_RANGE=2.0
INPUT_FILE="$WORK/lpa_normalized/LPA_ID_scale_zdepths.txt.gz"
OUTPUT_PREFIX="$WORK/lpa_neighbors/LPA_find_neighbor"

# Create output directory if it doesn't exist
OUTPUT_DIR=$(dirname "${OUTPUT_PREFIX}")
mkdir -p "${OUTPUT_DIR}"
cd "${OUTPUT_DIR}" || { echo "ERROR: Cannot change to output directory: ${OUTPUT_DIR}"; exit 1; }

# Print job information
echo "Job started at: $(date)"
echo "Running on node: $(hostname)"
echo "Job ID: $SLURM_JOB_ID"
echo ""
echo "Parameters:"
echo "  Batch number: ${BATCH_NUM}"
echo "  Total batches: ${TOTAL_BATCHES}"
echo "  Max z-range: ${MAX_Z_RANGE}"
echo "  Input file: ${INPUT_FILE}"
echo "  Output prefix: ${OUTPUT_PREFIX}"
echo ""

# Check if input file exists
if [[ ! -f "${INPUT_FILE}" ]]; then
    echo "ERROR: Input file does not exist: ${INPUT_FILE}"
    exit 1
fi

FIND_NEIGHBORS_PY="$LPA/playground/cram_crai/step5_find_neighbors/find_neighbors.py"
if [[ ! -f "${FIND_NEIGHBORS_PY}" ]]; then
    echo "ERROR: Python script not found: ${FIND_NEIGHBORS_PY}"
    exit 1
fi

INSPECT_FILE="$LPA/playground/cram_crai/step5_find_neighbors/inspect_file_format.py"
if [[ ! -f "${INSPECT_FILE}" ]]; then
    echo "ERROR: Inspect file script not found: ${INSPECT_FILE}"
    exit 1
fi
# Optional: Inspect input file format
echo "Inspecting input file format..."
python "${INSPECT_FILE}" "${INPUT_FILE}" > "inspect_input_file_log.txt"

# Run the Python script
echo "Starting find_neighbors analysis..."
python "${FIND_NEIGHBORS_PY}" \
    "${BATCH_NUM}" \
    "${TOTAL_BATCHES}" \
    "${MAX_Z_RANGE}" \
    "${INPUT_FILE}" \
    "${OUTPUT_PREFIX}" \
    > "running_log_find_neighbor_LPA_zmax=${MAX_Z_RANGE}.txt" 2>&1

# Check exit status
if [[ $? -eq 0 ]]; then
    echo "find_neighbors completed successfully!"
    echo "Output files:"
    ls -la "${OUTPUT_PREFIX}"*.txt.gz
else
    echo "ERROR: find_neighbors failed with exit code $?"
    exit 1
fi

echo ""
echo "Job completed at: $(date)"

# Optional: Print some statistics about the output
OUTPUT_FILE="${OUTPUT_PREFIX}.zMax${MAX_Z_RANGE}.txt.gz"
if [[ -f "${OUTPUT_FILE}" ]]; then
    echo ""
    echo "Output file statistics:"
    echo "  File: ${OUTPUT_FILE}"
    echo "  Size: $(ls -lh "${OUTPUT_FILE}" | awk '{print $5}')"
    echo "  Lines: $(zcat "${OUTPUT_FILE}" | wc -l)"
fi