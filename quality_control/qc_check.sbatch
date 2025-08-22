#!/bin/bash
#SBATCH --job-name=qc_check
#SBATCH --output=slurm/qc_check.out
#SBATCH --error=slurm/qc_check.err
#SBATCH --time=64:00:00
#SBATCH --partition=general
#SBATCH -n 1
#SBATCH --cpus-per-task=16
#SBATCH --mem=300g
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ztcaterer@colorado.edu

# Exit on any error, unset variable, or pipeline failure
set -euo pipefail

# Load necessary modules
module purge
module load samtools/1.21

# Define paths and filenames
VCF_SOURCE=$SOL/sequenced/cardiac_cohorts_SOL_dragen.chr6.vcf.gz
WORKDIR=$WORK/qc_check
OUTDIR=$LPA/quality_control/qc_check_output
VCF_NAME=$(basename "$VCF_SOURCE")
CLEANED_NAME=cardiac_cohorts_SOL_dragen_cleaned.chr6.vcf.gz

# Setup working directory
mkdir -p "$WORKDIR"
mkdir -p "$OUTDIR"
cd "$WORKDIR"

# Copy and confirm full file copy
echo "Copying VCF source to working directory..."
cp "$VCF_SOURCE" .
if [ ! -s "$VCF_NAME" ]; then
    echo "Copy failed or file is empty. Exiting."
    exit 1
fi

echo "Original VCF file size:"
ls -lh "$VCF_NAME"

# Optional: Count lines before cleaning
echo "Counting original lines..."
zcat "$VCF_NAME" | wc -l > linecount_before.txt

# Remove problematic line and compress safely
echo "Removing line 2336741 and compressing with bgzip..."
zcat "$VCF_NAME" | sed '2336741d' | bgzip -@ 16 -c > "$CLEANED_NAME"

# Check cleaned file integrity
echo "Verifying cleaned VCF..."
if ! bcftools view "$CLEANED_NAME" -Ov -o /dev/null; then
    echo "Cleaned VCF is malformed. Exiting."
    exit 1
fi

# Index the cleaned VCF
echo "Indexing cleaned VCF..."
tabix -p vcf "$CLEANED_NAME"

# Count lines after cleaning
zcat "$CLEANED_NAME" | wc -l > linecount_after.txt

# Move cleaned file and line counts to output
mv "$CLEANED_NAME" "$CLEANED_NAME.tbi" linecount_before.txt linecount_after.txt "$OUTDIR"

echo "✅ QC check completed successfully!"

# Cleanup
cd ..
rm -rf "$WORKDIR"
echo "Temporary files cleaned up."
echo "Output files moved to $OUTDIR"
echo "Job finished at $(date)"
