#!/bin/bash
#SBATCH --job-name=eagle_chr6_plink_LPA
#SBATCH --output=slurm/eagle_chr6_plink_LPA.out
#SBATCH --error=slurm/eagle_chr6_plink_LPA.err
#SBATCH --time=11-00:00:00
#SBATCH --partition=general
#SBATCH -n 1
#SBATCH --cpus-per-task=16
#SBATCH --mem=320g
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ztcaterer@colorado.edu

set -euo pipefail
echo "Job started on $(hostname) at $(date)"

# Load modules
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
WORKDIR=$WORK/eagle_chr6_plink_LPA
OUTDIR=$LPA/eagle_output_files/lpa_gene

mkdir -p "$WORKDIR" "$OUTDIR"
cd "$WORKDIR"

# Copy PLINK and map files
cp "${PLINK_ORIG_PREFIX}".bed .
cp "${PLINK_ORIG_PREFIX}".bim .
cp "${PLINK_ORIG_PREFIX}".fam .
cp "$GEN_MAP" .

# Filter region of interest (chr6:110Mb–210Mb) and assign SNP IDs
echo "🔍 Filtering LPA region and assigning unique SNP IDs..."
plink \
  --bfile cardiac_chr6 \
  --chr 6 \
  --from-bp 110000000 \
  --to-bp 210000000 \
  --set-missing-var-ids @:# \
  --make-bed \
  --out cardiac_chr6_LPAregion

# Run Eagle
echo "🚀 Running Eagle on LPA region..."
"$EAGLE_EXE" \
  --numThreads 16 \
  --bfile=cardiac_chr6_LPAregion \
  --chrom=6 \
  --geneticMapFile=$(basename "$GEN_MAP") \
  --outPrefix=phased_chr6_LPA_plink

# Move outputs
mv phased_chr6_LPA_plink.* "$OUTDIR"/

echo "✅ Eagle LPA phasing (PLINK) completed at $(date)"

# Cleanup
cd ..
rm -rf "$WORKDIR"
echo "✅ Cleanup completed at $(date)"
