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

GROUP="${FILEDIR}/sample_ids_split_${i}"
VCF_OUT="${FILEDIR}/subset_${i}.vcf.gz"

echo "🧬 Running gnomix on ${VCF_OUT}"
python $GNOMIX \
    "$VCF_OUT" \
    "$OUTDIR" \
    6 \
    True \
    # Genetic Map Files .tsv (w/ 3 columns chromosome #, snp physical location, and SNP genetic position) \
    # reference files (is a .vcf or .vcf.gz file containing reference haplotypes in any order) \
    # sample map file matching reference samples to their respectice reference populations

