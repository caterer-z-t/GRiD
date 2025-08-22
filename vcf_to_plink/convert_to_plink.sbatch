#!/bin/bash
#SBATCH --job-name=plink_chr6
#SBATCH --output=slurm/plink_chr6.out
#SBATCH --error=slurm/plink_chr6.err
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
WORKDIR=$WORK/plink_conv
OUTDIR=$HOME/LPA/vcf_to_plink/plink_output
VCF_SOURCE=$SOL/sequenced/cardiac_cohorts_SOL_dragen.chr6.vcf.gz
VCF_NAME=$(basename "$VCF_SOURCE")
CLEANED_VCF="cleaned_${VCF_NAME}"

# Setup work environment
mkdir -p "$WORKDIR"
mkdir -p "$OUTDIR"
cd "$WORKDIR"
echo "Working in $WORKDIR"

# Copy VCF
cp "$VCF_SOURCE" .

# Clean VCF by removing the malformed line
echo "Cleaning VCF (removing line 2336741)..."
zcat "$VCF_NAME" | sed '2336741d' | bgzip -c > "$CLEANED_VCF"
tabix -p vcf "$CLEANED_VCF"

# Run PLINK
echo "Running plink2 on cleaned VCF..."
plink2 \
  --vcf "$CLEANED_VCF" \
  --make-bed \
  --out cardiac_chr6 \
  --threads 16

# Move results
mv cardiac_chr6* "$OUTDIR"

# Cleanup
cd ..
rm -rf "$WORKDIR"

echo "✅ Chromosome 6 PLINK conversion complete."
