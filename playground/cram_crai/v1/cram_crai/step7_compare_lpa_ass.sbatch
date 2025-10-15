#!/bin/bash
#SBATCH --job-name=plot_vntr_lpa
#SBATCH --output=slurm/%x.out
#SBATCH --error=slurm/%x.err
#SBATCH --time=01:00:00
#SBATCH --partition=general
#SBATCH -n 1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ztcaterer@colorado.edu

source "$GRID/playground/cram_crai/common.sh"

usage() {
    echo "Usage: $0 [-p PYTHON_SCRIPT] [-d DIPLOID_FILE] [-f PHENO_FILE] [-o OUT_PREFIX] [-O OUT_DIR] [-h]"
    echo ""
    echo "  -p   Python script to run"
    echo "  -d   Diploid CN estimate file (two columns: nwd_id, value)"
    echo "  -f   Phenotype CSV file (contains nwd_id and lpa_nmol)"
    echo "  -o   Output prefix for plots and results"
    echo "  -O  Output directory"
    echo "  -h   Show this help message"
    exit 1
}

# Parse command line arguments
while getopts "p:d:f:o:O:h" opt; do
    case $opt in
        p) PYTHON_SCRIPT="$OPTARG" ;;
        d) DIPLOID_FILE="$OPTARG" ;;
        f) PHENO_FILE="$OPTARG" ;;
        o) OUT_PREFIX="$OPTARG" ;;
        O) OUT_DIR="$OPTARG" ;;
        h) usage ;;
        *) usage ;;
    esac
done

require_args PYTHON_SCRIPT DIPLOID_FILE PHENO_FILE OUT_PREFIX OUT_DIR

export PYTHON_SCRIPT DIPLOID_FILE PHENO_FILE OUT_PREFIX OUT_DIR

# -------------------------------
# Modules / Environment
# -------------------------------
module purge
module load anaconda
source $(conda info --base)/etc/profile.d/conda.sh
eval "$(conda shell.bash hook)"
conda activate lpa_vntr

# -------------------------------
# Validate Inputs
# -------------------------------
require_file "$PYTHON_SCRIPT"
require_file "$DIPLOID_FILE"
require_file "$PHENO_FILE"

mkdir -p $OUT_DIR
cd $OUT_DIR

# -------------------------------
# Run Plotting Script
# -------------------------------
python "$PYTHON_SCRIPT" \
  --diploid_file "$DIPLOID_FILE" \
  --pheno_file "$PHENO_FILE" \
  --output_prefix "$OUT_PREFIX" \
  > "${OUT_PREFIX}_run.log" 2>&1

echo "Plotting complete. Outputs saved with prefix: $OUT_PREFIX"
