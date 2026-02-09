#!/bin/bash
#SBATCH --job-name=IBS_IBD
#SBATCH --output=slurm/%x.out
#SBATCH --error=slurm/%x.err
#SBATCH --time=2-00:00:00
#SBATCH --mem=320G
#SBATCH --partition=general
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ztcaterer@colorado.edu

set -euo pipefail

####################
# Modules
####################
module purge
module load plink/1.90b6.21
module load gcc boost zlib

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

####################
# Directories
####################
WORKDIR=$WORK/IBS_IBD1_IBD2
mkdir -p "$WORKDIR"/{plink,cpp,logs}
cd "$WORKDIR"

####################
# Inputs
####################
GNOMIX_VCF="$WORK/gnomix/files/phased_samples_gnomix.vcf.gz"
PLINK_PREFIX=plink/phased_samples_gnomix

####################
# Step 1: VCF → BED
####################
plink \
  --vcf "$GNOMIX_VCF" \
  --make-bed \
  --out "$PLINK_PREFIX"

####################
# Step 2: PLINK IBS + IBD
####################
plink \
  --bfile "$PLINK_PREFIX" \
  --genome full \
  --threads $SLURM_CPUS_PER_TASK \
  --out "$PLINK_PREFIX"

####################
# Step 3: Extract IBD2 and IBD1
####################
GENOME_FILE="$PLINK_PREFIX.genome"

# IBD2 (Z2 ≥ 0.75)
awk 'NR==1 || $10 >= 0.75' "$GENOME_FILE" \
  > plink/IBD2_pairs.txt

# IBD1 (0.25 ≤ Z1 ≤ 0.75)
awk 'NR==1 || ($9 >= 0.25 && $9 <= 0.75)' "$GENOME_FILE" \
  > plink/IBD1_pairs.txt

####################
# Step 4: IBS via PLINK
####################
# IBS0 is column 8
awk 'NR==1 || $8 > 0' "$GENOME_FILE" \
  > plink/IBS_plink.txt

####################
# Step 5: IBS via C++ code
####################
VCF_CPP="$GNOMIX_VCF"
OUT_CPP=cpp/IBS_custom.txt.gz

VNTR_START=123456
VNTR_END=123999

bin/compute_IBS \
  "$VCF_CPP" \
  "$OUT_CPP" \
  "$VNTR_START" \
  "$VNTR_END" \
  0.5 \
  0.5 \
  $SLURM_CPUS_PER_TASK

echo "✅ PLINK IBS/IBD + custom IBS complete"
