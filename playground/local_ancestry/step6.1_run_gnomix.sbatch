#!/bin/bash
#SBATCH --job-name=gnomix
#SBATCH --output=slurm/gnomix_%a.out
#SBATCH --error=slurm/gnomix_%a.err
#SBATCH --time=11-00:00:00
#SBATCH --mem=300G
#SBATCH -n 1
#SBATCH --cpus-per-task=64
#SBATCH --array=0-10   # <--- Change this depending on how many array tasks you want
#SBATCH --mail-type=ALL 
#SBATCH --mail-user=ztcaterer@colorado.edu

set -euo pipefail

module purge
module load anaconda 
source $(conda info --base)/etc/profile.d/conda.sh
conda activate gnomix
module load bcftools

export NUMEXPR_MAX_THREADS=$SLURM_CPUS_PER_TASK

WORKDIR=$WORK/gnomix
FILEDIR=$WORKDIR/files
OUTDIR=$WORKDIR/output
mkdir -p $OUTDIR

cd $WORKDIR

GNOMIX=$WORKDIR/gnomix.py
IN_VCF="${FILEDIR}/phased_merged_SOL_ref.vcf.gz"
[ ! -f "$IN_VCF" ] && IN_VCF="${IN_VCF%.gz}"

# total number of chunks
TOTAL_SPLITS=100

# number of array tasks (auto-detected from sbatch --array)
NUM_TASKS=$SLURM_ARRAY_TASK_COUNT

# how many splits each task should process
CHUNK_SIZE=$(( (TOTAL_SPLITS + NUM_TASKS - 1) / NUM_TASKS ))

# range for this array task
START=$(( SLURM_ARRAY_TASK_ID * CHUNK_SIZE ))
END=$(( START + CHUNK_SIZE - 1 ))

# cap at TOTAL_SPLITS-1
if [ $END -ge $TOTAL_SPLITS ]; then
    END=$(( TOTAL_SPLITS - 1 ))
fi

echo "Task ${SLURM_ARRAY_TASK_ID} handling splits $START to $END"

# check for pretrained models dir (same as before)
if [ ! -d "$WORKDIR/pretrained_gnomix_models" ]; then
    echo "Pretrained models not found, downloading..."
    if ! command -v gdown &> /dev/null; then
        pip install --user gdown
        export PATH=$PATH:$HOME/.local/bin
    fi
    gdown https://drive.google.com/uc?id=1Q0zg9zqTaZUvt42uxzE_0gcfnFvjwi33 -O pretrained_gnomix_models.tar.gz
    tar -xvf pretrained_gnomix_models.tar.gz
    rm pretrained_gnomix_models.tar.gz
fi

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
