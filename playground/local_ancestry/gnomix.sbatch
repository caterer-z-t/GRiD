#!/bin/bash
#SBATCH --job-name=gnomix
#SBATCH --output=slurm/%x.out
#SBATCH --error=slurm/%x.err
#SBATCH --time=11-00:00:00
#SBATCH --mem=300G
#SBATCH -n 1
#SBATCH --cpus-per-task=64
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ztcaterer@colorado.edu

# set -euo pipefail

############################
# Environment
############################
module purge
module load anaconda
module load bcftools

source $(conda info --base)/etc/profile.d/conda.sh
eval "$(conda shell.bash hook)"
conda activate gnomix

export NUMEXPR_MAX_THREADS=64
export OMP_NUM_THREADS=64
export MKL_NUM_THREADS=64

############################
# Directories
############################
WORKDIR=$WORK/gnomix_v2
FILEDIR=$WORKDIR/files
OUTDIR=$WORKDIR/output

mkdir -p "$WORKDIR" "$FILEDIR" "$OUTDIR"
cd "$WORKDIR"

############################
# Gnomix
############################
if [ ! -d "$WORKDIR/gnomix" ]; then
    cp -r "$SOFTWARE/gnomix" "$WORKDIR"
fi
GNOMIX="$WORKDIR/gnomix/gnomix.py"

############################
# Input VCF
############################
IN_VCF="${SOL}/sequenced/cardiac_cohorts_SOL_dragen.chr6.vcf.gz"

############################
# Detect chromosome naming
############################
echo "🔍 Detecting chromosome naming..."
if bcftools view -h "$IN_VCF" | grep -q "##contig=<ID=chr6"; then
    CHR="chr6"
elif bcftools view -h "$IN_VCF" | grep -q "##contig=<ID=6"; then
    CHR="6"
else
    echo "❌ ERROR: chr6 not found in VCF header"
    bcftools view -h "$IN_VCF" | grep "^##contig" | head
    CHR="chr6"
fi

echo "✅ Using chromosome: $CHR"

############################
# Subset region (LPA)
############################
SUB_VCF_RAW="${FILEDIR}/cardiac_cohorts_SOL_dragen.${CHR}.subset.raw.vcf.gz"
SUB_VCF="${FILEDIR}/cardiac_cohorts_SOL_dragen.${CHR}.subset.vcf.gz"

# Subset region
bcftools view \
    -r ${CHR}:135000000-185000000 \
    -Oz \
    -o "$SUB_VCF_RAW" \
    "$IN_VCF" \
    --threads 64

# Normalize and clean for Eagle
bcftools norm \
    -m -any \
    "$SUB_VCF_RAW" | \
bcftools view \
    -m2 -M2 \
    -v snps \
    --min-ac 1 \
    -i 'F_MISSING < 0.1' | \
bcftools +setGT -- -t q -n . -i 'GT!="0/0" & GT!="0/1" & GT!="1/0" & GT!="1/1" & GT!="./." & GT!="0|0" & GT!="0|1" & GT!="1|0" & GT!="1|1" & GT!=".|."' | \
bcftools view \
    -Oz \
    -o "$SUB_VCF" \
    --threads 64

bcftools index "$SUB_VCF"

############################
# Check variant count
############################
NVAR=$(bcftools view -H "$SUB_VCF" | wc -l)

if [ "$NVAR" -eq 0 ]; then
    echo "❌ ERROR: Subset VCF contains ZERO variants"
    exit 1
fi

echo "✅ Subset contains $NVAR variants ready for phasing"

############################
# Eagle phasing
############################
EAGLE="${SOFTWARE}/Eagle_v2.4.1/eagle"
EAGLE_DIR="${SOFTWARE}/Eagle_v2.4.1"

PHASED_PREFIX="${FILEDIR}/cardiac_cohorts_SOL_dragen.${CHR}.subset.phased"
PHASED_VCF="${PHASED_PREFIX}.vcf.gz"

echo "🧬 Running Eagle phasing..."

$EAGLE \
    --vcf "$SUB_VCF" \
    --geneticMapFile "${EAGLE_DIR}/tables/genetic_map_hg38_withX.txt.gz" \
    --outPrefix "$PHASED_PREFIX" \
    --numThreads 64

############################
# Check phased VCF
############################
if [ ! -f "$PHASED_VCF" ]; then
    echo "❌ ERROR: Eagle did not produce phased VCF"
    exit 1
fi

############################
# Pretrained Gnomix models
############################
cd "$WORKDIR"

if [ ! -d "$WORKDIR/pretrained_gnomix_models" ]; then
    echo "⬇️ Downloading pretrained Gnomix models..."
    pip install --user gdown
    export PATH=$PATH:$HOME/.local/bin

    gdown https://drive.google.com/uc?id=1Q0zg9zqTaZUvt42uxzE_0gcfnFvjwi33 \
        -O pretrained_gnomix_models.tar.gz

    tar -xvf pretrained_gnomix_models.tar.gz
    rm pretrained_gnomix_models.tar.gz
fi

############################
# Run Gnomix
############################
echo "🚀 Running Gnomix..."

python "$GNOMIX" \
    "$PHASED_VCF" \
    "$OUTDIR" \
    "$CHR" \
    true \
    "$WORKDIR/pretrained_gnomix_models/${CHR}/model_chm_6.pkl"

echo "✅ Gnomix complete"
