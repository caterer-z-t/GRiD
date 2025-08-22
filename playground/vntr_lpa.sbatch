#!/bin/bash
#SBATCH --job-name=vntr_lpa
#SBATCH --output=slurm/vntr_lpa.out
#SBATCH --error=slurm/vntr_lpa.err
#SBATCH --time=11-00:00:00
#SBATCH --partition=general
#SBATCH -n 1
#SBATCH --cpus-per-task=16
#SBATCH --mem=300g
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ztcaterer@colorado.edu

# Exit on error and undefined variables
set -euo pipefail

echo "Job started on $(hostname) at $(date)"

# Load necessary modules
module purge

# Define paths
CLEANED_VCF=$HOME/epi_caterer/lpa/playground/qc_check_output/cardiac_cohorts_SOL_dragen_cleaned.chr6.vcf.gz
WORKDIR=$WORK/vntr_lpa
OUTDIR=$HOME/epi_caterer/lpa/playground/vntr_lpa_output

# Prepare working directory
mkdir -p "$WORKDIR" "$OUTDIR"
cd "$WORKDIR"

# Copy required files
cp "$CLEANED_VCF" .
cp "${CLEANED_VCF}.tbi" .

