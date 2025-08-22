#!/bin/bash
#SBATCH --job-name=local_ancestry
#SBATCH --output=slurm/local_ancestry.out
#SBATCH --error=slurm/local_ancestry.err
#SBATCH --time=11-00:00:00
#SBATCH --mem=500G
#SBATCH -n 1
#SBATCH --cpus-per-task=2
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ztcaterer@colorado.edu

set -euo pipefail
echo "Job started on $(hostname) at $(date)"

# Load modules
module purge
module load anaconda
module load bcftools

# Initialize Conda (very important in batch scripts)
source $(conda info --base)/etc/profile.d/conda.sh
eval "$(conda shell.bash hook)"

source /nas/longleaf/apps/anaconda/4.3.0/anaconda/etc/profile.d/conda.sh
conda activate gnomix

# Work directory
WORKDIR=$WORK/gnomix
FILEDIR=$WORKDIR/files
OUTDIR=$WORK/gnomix/output
mkdir -p $WORKDIR $FILEDIR $OUTDIR
cd $WORKDIR

CH6_PREFIX=phased_merged_SOL_ref
LPA_PREFIX=phased_merged_SOL_ref_LPA

# Define source and destination directories
SOURCE_DIR=$WORK/phasing/

# Array of prefixes
prefixes=("${CH6_PREFIX}" "${LPA_PREFIX}")

for prefix in "${prefixes[@]}"; do
    if [ -f "${SOURCE_DIR}/${prefix}.vcf" ]; then
        echo "📄 Copying ${prefix}.vcf"
        cp "${SOURCE_DIR}/${prefix}.vcf" "${FILEDIR}/"
    elif [ -f "${SOURCE_DIR}/${prefix}.vcf.gz" ]; then
        echo "📄 Copying ${prefix}.vcf.gz"
        cp "${SOURCE_DIR}/${prefix}.vcf.gz" "${FILEDIR}/"
    else
        echo "⚠️  No VCF file found for ${prefix}"
    fi
done

# Copy gnomix into the working directory
if [ ! -f "${WORKDIR}/gnomix.py" ]; then
    echo "📁 Copying gnomix directory"
    cp -r ${SOFTWARE}/gnomix/* ${WORKDIR}/
else
    echo "📁 gnomix directory already exists in ${WORKDIR}"
fi

# Extract and split sample IDs
VCF_FILE="${FILEDIR}/${CH6_PREFIX}.vcf.gz"
[ ! -f "$VCF_FILE" ] && VCF_FILE="${FILEDIR}/${CH6_PREFIX}.vcf"

bcftools query -l "$VCF_FILE" > "${FILEDIR}/sample_ids.txt"

split -n l/100 -d -a 2 "${FILEDIR}/sample_ids.txt" "${FILEDIR}/sample_ids_split_"

# Path to gnomix.py
GNOMIX=gnomix.py

# python $GNOMIX <query_file> <output_folder> <chr_nr> <phase> <path_to_model> 

# advanced gnomix options
# $ python3 gnomix.py \
#    <query_file> \
#    <output_folder> \
#    <chr_nr> \
#    <phase> \
#    <genetic_map_file> \
#    <reference_file> \
#    <sample_map_file> \
#    <config_file> 


# Parse optional region flag
REGION=${1:-"all"}
LPA_FLAG=false

if [[ "$REGION" =~ ^(-y|-yes|-Y|-YES|-lpa|-LPA)$ ]]; then
    LPA_FLAG=true
fi

# Gnomix execution
for group in "${FILEDIR}/sample_ids_split_"*; do
    BASENAME=$(basename "$group")
    VCF_OUT="${FILEDIR}/subset_${BASENAME}.vcf.gz"

    if $LPA_FLAG; then
        IN_VCF="${FILEDIR}/${LPA_PREFIX}.vcf.gz"
    else
        IN_VCF="${FILEDIR}/${CH6_PREFIX}.vcf.gz"
    fi

    [ ! -f "$IN_VCF" ] && IN_VCF="${IN_VCF%.gz}"  # fallback to uncompressed

    echo "📤 Subsetting VCF for ${BASENAME}"
    bcftools view -S "$group" -Oz -o "$VCF_OUT" "$IN_VCF"
    bcftools index "$VCF_OUT"

    echo "🧬 Running gnomix on ${VCF_OUT}"
    python $GNOMIX \
        "$VCF_OUT" \
        "$OUTDIR" \
        6 \
        True \
        "$WORKDIR/pretrained_gnomix_models/chr6/model_chm_6.pkl"

    # rename output files
    mv "$OUTDIR/query_file_phased.vcf" "${OUTDIR}/gnomix_output_${BASENAME}.vcf"
    mv "$OUTDIR/query_results.fb" "${OUTDIR}/gnomix_results_${BASENAME}.fb"
    mv "$OUTDIR/query_results.msp" "${OUTDIR}/gnomix_msp_${BASENAME}.msp"
    
done