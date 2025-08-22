#!/bin/bash
#SBATCH --job-name=vcf2plink
#SBATCH --output=slurm/vcf2plink.out
#SBATCH --error=slurm/vcf2plink.err
#SBATCH --time=2-00:00:00
#SBATCH --partition=general
#SBATCH -n 1
#SBATCH --cpus-per-task=5
#SBATCH --mem=500G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ztcaterer@colorado.edu

# Load required modules
module purge
module load samtools/1.20
module load plink/1.90b6.21

# Set working directory
WORKDIR=$WORK/vcf2plink
mkdir -p "$WORKDIR"
cd "$WORKDIR" || exit 1

# Define paths to input VCFs
SOL_CH6="$SOL/sequenced/cardiac_cohorts_SOL_dragen.chr6.vcf.gz"
REFERENCE_CH6="$LPA/reference_panel/gnomad.genomes.v3.1.2.hgdp_tgp.chr6.vcf.bgz"

# Get basenames
SOL_BASENAME=$(basename "$SOL_CH6")
REFERENCE_BASENAME=$(basename "$REFERENCE_CH6")

# Copy input files if not present
[ ! -f "$SOL_BASENAME" ] && cp "$SOL_CH6" .
[ ! -f "$REFERENCE_BASENAME" ] && cp "$REFERENCE_CH6" .

CLEANED_SOL_VCF=$LPA/quality_control/qc_check_output/cardiac_cohorts_SOL_dragen_cleaned.chr6.vcf.gz
[ ! -f "$CLEANED_SOL_VCF" ] && cp "$CLEANED_SOL_VCF" .

# Unzip cleaned SOL if needed
CLEANED_SOL_VCF_UNZIPPED="${CLEANED_SOL_VCF%.gz}"
[ ! -f "$CLEANED_SOL_VCF_UNZIPPED" ] && gunzip -f "$CLEANED_SOL_VCF"

# Decompress reference panel VCF if needed
REFERENCE_VCF="${REFERENCE_BASENAME%.bgz}"
[ ! -f "$REFERENCE_VCF" ] && bgzip -d -f "$REFERENCE_BASENAME"

# Confirm input files
echo "Input VCFs ready:"
ls -lh "$CLEANED_SOL_VCF_UNZIPPED" "$REFERENCE_VCF"

# Convert SOL VCF to PLINK
if [ ! -f SOL_ch6.bed ]; then
    echo "Converting cleaned SOL VCF to PLINK format..."
    plink \
        --vcf "$CLEANED_SOL_VCF_UNZIPPED" \
        --double-id \
        --make-bed \
        --out SOL_ch6 \
        --threads 5
else
    echo "SOL_ch6.bed already exists. Skipping PLINK conversion for SOL."
fi

# Convert reference panel to PLINK
if [ ! -f reference_panel.bed ]; then
    echo "Converting reference panel VCF to PLINK format..."
    plink \
        --vcf "$REFERENCE_VCF" \
        --double-id \
        --make-bed \
        --out reference_panel \
        --threads 5
else
    echo "reference_panel.bed already exists. Skipping PLINK conversion for reference."
fi

echo "All steps completed successfully."
