#!/bin/bash
#SBATCH --job-name=plink_chr6_cleaned
#SBATCH --output=slurm/plink_chr6_cleaned.out
#SBATCH --error=slurm/plink_chr6_cleaned.err
#SBATCH --time=11-00:00:00
#SBATCH --partition=general
#SBATCH -n 1
#SBATCH --cpus-per-task=16
#SBATCH --mem=300g
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ztcaterer@colorado.edu

# Exit on error
set -euo pipefail

# Load modules
module purge
module load plink/2.00a
module load samtools/1.21

# Set paths and filenames
WORKDIR=$WORK/plink_conv_cleaned
OUTDIR=$LPA/vcf_to_plink/plink_output_cleaned
CLEANED_VCF=$LPA/quality_control/qc_check_output/cardiac_cohorts_SOL_dragen_cleaned.chr6.vcf.gz

# Setup work environment
mkdir -p "$WORKDIR"
cd "$WORKDIR"
echo "Working in $WORKDIR"

# Copy cleaned VCF and index
cp "$CLEANED_VCF" .
cp "${CLEANED_VCF}.tbi" .

# Run PLINK
echo "Running plink2 on cleaned VCF..."
plink2 \
  --vcf $(basename "$CLEANED_VCF") \
  --make-bed \
  --out cardiac_chr6 \
  --threads 16

# Move results
mkdir -p "$OUTDIR"
mv cardiac_chr6* "$OUTDIR"

echo "✅ PLINK conversion completed at $(date)"

# Clean up
echo "Cleaning up..."
cd ..
rm -rf "$WORKDIR"
echo "✅ Cleanup completed at $(date)"