#!/bin/bash
#SBATCH --job-name=phase
#SBATCH --output=slurm/%x.out
#SBATCH --error=slurm/%x.err
#SBATCH --time=11-00:00:00
#SBATCH --mem=320G
#SBATCH --partition=general
#SBATCH -n 1
#SBATCH --cpus-per-task=10
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ztcaterer@colorado.edu

set -euo pipefail
echo "Job started on $(hostname) at $(date)"

##################################
# Load modules
##################################
module purge
module load anaconda
module load plink/1.90b6.21

# Initialize conda ONCE
source $(conda info --base)/etc/profile.d/conda.sh
conda activate lpa

##################################
# Paths
##################################
EAGLE_EXE="$SOFTWARE/Eagle_v2.4.1/eagle"
GEN_MAP="$SOFTWARE/Eagle_v2.4.1/tables/genetic_map_1cMperMb.txt"

PLINK_PREFIX="$WORK/merge_files/merged_SOL_ref"

WORKDIR="$WORK/phasing"
mkdir -p "$WORKDIR"
cd "$WORKDIR"

##################################
# Copy PLINK files
##################################
for ext in bed bim fam; do
    cp "${PLINK_PREFIX}.${ext}" .
done

cp "$GEN_MAP" .

PREFIX=$(basename "$PLINK_PREFIX")


##################################
# 2. Run Eagle
##################################
echo "🚀 Running Eagle phasing..."

plink --bfile "$PREFIX" --recode vcf --out "$PREFIX"

"$EAGLE_EXE" \
  --numThreads 10 \
  --vcf "$PREFIX.vcf" \
  --geneticMapFile "$(basename "$GEN_MAP")" \
  --chrom 6 \
  --vcfOutFormat v \
  --outPrefix phased_merged_SOL_ref

echo "✅ Phasing completed at $(date)"
