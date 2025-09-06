#!/bin/bash
#SBATCH --job-name=gnomix
#SBATCH --output=slurm/gnomix_%A_%a.out
#SBATCH --error=slurm/gnomix_%A_%a.err
#SBATCH --time=5-00:00:00
#SBATCH --mem=200G
#SBATCH -n 1
#SBATCH --cpus-per-task=2
#SBATCH --array=0-3   # 4 jobs in parallel
#SBATCH --mail-type=ALL 
#SBATCH --mail-user=ztcaterer@colorado.edu

set -euo pipefail

module purge
module load anaconda 
source $(conda info --base)/etc/profile.d/conda.sh
conda activate gnomix
module load bcftools

# list the given libraries in the conda environment
conda list

WORKDIR=$WORK/gnomix
FILEDIR=$WORKDIR/files
OUTDIR=$WORKDIR/output
mkdir -p $OUTDIR

cd $WORKDIR

GNOMIX=$WORKDIR/gnomix.py
IN_VCF="${FILEDIR}/phased_merged_SOL_ref.vcf.gz"
[ ! -f "$IN_VCF" ] && IN_VCF="${IN_VCF%.gz}"

# Which slice of splits does this array task handle?
START=$(( SLURM_ARRAY_TASK_ID * 25 ))
END=$(( START + 24 ))

# check for pretrained models dir
if [ ! -d "$WORKDIR/pretrained_gnomix_models" ]; then
    echo "Pretrained models not found, downloading..."
    
    # Install gdown if not already available
    if ! command -v gdown &> /dev/null; then
        echo "gdown not found, installing via pip..."
        pip install --user gdown
        export PATH=$PATH:$HOME/.local/bin
    fi

    # Download pretrained models from Google Drive
    gdown https://drive.google.com/uc?id=1Q0zg9zqTaZUvt42uxzE_0gcfnFvjwi33 -O pretrained_gnomix_models.tar.gz
    
    echo "Expanding pretrained models..."
    tar -xvf pretrained_gnomix_models.tar.gz
    rm pretrained_gnomix_models.tar.gz

    echo "Downloaded pretrained models to directory ./pretrained_gnomix_models"
else
    echo "Pretrained models found, skipping download."
fi

echo "Handling splits $START to $END"

for i in $(seq -w $START $END); do
    GROUP="${FILEDIR}/sample_ids_split_${i}"
    VCF_OUT="${FILEDIR}/subset_${i}.vcf.gz"

    echo "📤 Subsetting VCF for ${i}"
    bcftools view -S "$GROUP" -Oz -o "$VCF_OUT" "$IN_VCF"
    bcftools index "$VCF_OUT"

    echo "🧬 Running gnomix on ${VCF_OUT}"
    python $GNOMIX \
        "$VCF_OUT" \
        "$OUTDIR" \
        6 \
        True \
        "$WORKDIR/pretrained_gnomix_models/chr6/model_chm_6.pkl"

    mv "$OUTDIR/query_file_phased.vcf" "${OUTDIR}/gnomix_output_${i}.vcf"
    mv "$OUTDIR/query_results.fb" "${OUTDIR}/gnomix_results_${i}.fb"
    mv "$OUTDIR/query_results.msp" "${OUTDIR}/gnomix_msp_${i}.msp"
done
