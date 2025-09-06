#!/bin/bash
#SBATCH --job-name=split_samples
#SBATCH --output=slurm/split_samples.out
#SBATCH --error=slurm/split_samples.err
#SBATCH --time=02:00:00
#SBATCH --mem=50G
#SBATCH -n 1
#SBATCH --cpus-per-task=2
#SBATCH --mail-type=ALL 
#SBATCH --mail-user=ztcaterer@colorado.edu

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
        echo "⚠️ No VCF file found for ${prefix}"
    fi
done

# Copy gnomix into the working directory
if [ ! -f "${WORKDIR}/gnomix.py" ]; then
    echo "📁 Copying gnomix directory"
    cp -r ${SOFTWARE}/gnomix/* ${WORKDIR}/
else
    echo "📁 gnomix directory already exists in ${WORKDIR}"
fi

VCF_FILE="${FILEDIR}/phased_merged_SOL_ref.vcf.gz"
[ ! -f "$VCF_FILE" ] && VCF_FILE="${FILEDIR}/phased_merged_SOL_ref.vcf"

bcftools query -l "$VCF_FILE" > "${FILEDIR}/sample_ids.txt"

# Split into 100 files, sample_ids_split_00 .. 99
split -n l/100 -d -a 2 "${FILEDIR}/sample_ids.txt" "${FILEDIR}/sample_ids_split_"

