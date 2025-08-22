#!/bin/bash
#SBATCH --job-name=phase
#SBATCH --output=slurm/phase.out
#SBATCH --error=slurm/phase.err
#SBATCH --time=11-00:00:00
#SBATCH --mem=320G
#SBATCH --partition=general
#SBATCH -n 1
#SBATCH --cpus-per-task=10
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ztcaterer@colorado.edu

set -euo pipefail
echo "Job started on $(hostname) at $(date)"

# Load modules
module purge
module load anaconda

# Initialize Conda (very important in batch scripts)
source $(conda info --base)/etc/profile.d/conda.sh
eval "$(conda shell.bash hook)"

source /nas/longleaf/apps/anaconda/4.3.0/anaconda/etc/profile.d/conda.sh
conda activate lpa

module load plink/1.90b6.21

# Define paths for software
EAGLE_EXE=$SOFTWARE/Eagle_v2.4.1/eagle
GEN_MAP=$SOFTWARE/Eagle_v2.4.1/tables/genetic_map_1cMperMb.txt

# Merged PLINK input file prefix
PLINK_PATH=$WORK/pca/merged_SOL_ref
PREFIX=$(basename "$PLINK_PATH")
LPA_PRUNED_PREFIX="${PREFIX}_LPA"


# Define directories
WORKDIR=$WORK/phasing
mkdir -p "$WORKDIR"
cd "$WORKDIR"

# Copy PLINK and map files if they don't already exist
for ext in bed bim fam; do
    if [ ! -f "$PREFIX.${ext}" ]; then
        cp "${PLINK_PATH}.${ext}" .
    fi
done

if [ ! -f "$(basename "$GEN_MAP")" ]; then
    cp "$GEN_MAP" .
fi

# Parse optional region flag
REGION=${1:-"all"}
LPA_FLAG=false

if [[ "$REGION" =~ ^(-y|-yes|-Y|-YES|-lpa|-LPA)$ ]]; then
  LPA_FLAG=true
fi

# Run PLINK filtering and Eagle depending on flag
if $LPA_FLAG; then
  echo "🔍 Filtering LPA region (chr6:110–210Mb) and assigning SNP IDs..."

  # Check if the merged PLINK file for LPA region already exists
  if [ ! -f "merged_SOL_ref_LPA.bed" ]; then
    echo "📂 Extracting LPA region from PLINK file..."
    plink \
      --bfile $PREFIX \
      --chr 6 \
      --from-bp 110000000 \
      --to-bp 210000000 \
      --set-missing-var-ids @:# \
      --make-bed \
      --out $LPA_PREFIX \
      --threads 16
  fi

  # Check if the VCF file already exists
  if [ ! -f "$LPA_PREFIX.vcf" ]; then
    echo "🧬 Converting PLINK to VCF..."
    plink \
        --bfile $LPA_PREFIX \
        --recode vcf-iid \
        --chr 6 \
        --out $LPA_PREFIX \
        --threads 16
    fi 

  # Run Eagle on the LPA region
  echo "🚀 Running Eagle on LPA region using VCF input..."
  "$EAGLE_EXE" \
    --numThreads 16 \
    --vcf "$LPA_PREFIX.vcf" \
    --geneticMapFile=$(basename "$GEN_MAP") \
    --chrom 6 \
    --vcfOutFormat v \
    --outPrefix=phased_merged_SOL_ref_LPA

  echo "✅ LPA phasing completed at $(date)"
else
  # Check if vcf files aleady exist
  if [ ! -f "$PREFIX.vcf" ]; then
    echo "🧬 Converting PLINK to VCF..."
    plink \
      --bfile $PREFIX \
      --recode vcf-iid \
      --chr 6 \
      --out "$PREFIX" \
      --threads 16
  fi

  # Run Eagle on the whole chromosome
  echo "🚀 Running Eagle on whole chromosome using VCF input..."
  "$EAGLE_EXE" \
    --numThreads 16 \
    --vcf "$PREFIX.vcf" \
    --geneticMapFile=$(basename "$GEN_MAP") \
    --chrom 6 \
    --vcfOutFormat v \
    --outPrefix=phased_merged_SOL_ref

  echo "✅ Whole chromosome phasing completed at $(date)"
fi
