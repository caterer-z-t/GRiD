#!/bin/bash
#SBATCH --job-name=eagle_chr6_plink
#SBATCH --output=slurm/eagle_chr6_plink.out
#SBATCH --error=slurm/eagle_chr6_plink.err
#SBATCH --time=11-00:00:00
#SBATCH --partition=general
#SBATCH -n 1
#SBATCH --cpus-per-task=16
#SBATCH --mem=600g
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ztcaterer@colorado.edu

# Exit on error and undefined variables
set -euo pipefail

echo "Job started on $(hostname) at $(date)"

# Load necessary modules
module purge
module load anaconda/4.3.0
module load plink/1.90b6.21
source /nas/longleaf/apps/anaconda/4.3.0/anaconda/etc/profile.d/conda.sh
conda activate lpa

# Define paths for software
EAGLE_EXE=$SOFTWARE/Eagle_v2.4.1/eagle
GEN_MAP=$SOFTWARE/Eagle_v2.4.1/tables/genetic_map_1cMperMb.txt

# Plink output prefix
PLINK_ORIG_PREFIX=$LPA/vcf_to_plink/plink_output_cleaned/cardiac_chr6

# Define directories
WORKDIR=$WORK/eagle_chr6_plink
OUTDIR=$LPA/eagle_output_files/ch6

# Prepare working directory
mkdir -p "$WORKDIR" "$OUTDIR"
cd "$WORKDIR"

# Copy PLINK and genetic map files
cp "${PLINK_ORIG_PREFIX}".bed .
cp "${PLINK_ORIG_PREFIX}".bim .
cp "${PLINK_ORIG_PREFIX}".fam .
cp "$GEN_MAP" .

# Fix variant IDs
echo "🔧 Assigning unique SNP IDs with PLINK..."
plink \
  --bfile cardiac_chr6 \
  --set-missing-var-ids @:# \
  --make-bed \
  --out cardiac_chr6_fixed_ids

# Run Eagle using PLINK input
echo "🚀 Running Eagle phasing..."
"$EAGLE_EXE" \
  --numThreads 16 \
  --bfile=cardiac_chr6_fixed_ids \
  --chrom=6 \
  --geneticMapFile=$(basename "$GEN_MAP") \
  --outPrefix=phased_chr6_plink

# Move output files to final output directory
mv phased_chr6_plink.* "$OUTDIR"/

echo "✅ Eagle phasing (PLINK) completed at $(date)"

# Cleanup
cd ..
rm -rf "$WORKDIR"
echo "✅ Cleanup completed at $(date)"
