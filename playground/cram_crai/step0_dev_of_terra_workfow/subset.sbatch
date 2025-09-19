#!/bin/bash
#SBATCH --job-name=subset
#SBATCH --output=slurm/%x.out
#SBATCH --error=slurm/%x.err
#SBATCH --time=2-00:00:00
#SBATCH --partition=general
#SBATCH -n 1
#SBATCH --cpus-per-task=5
#SBATCH --mem=5G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ztcaterer@colorado.edu

# Load Modules
module purge
module load anaconda
source $(conda info --base)/etc/profile.d/conda.sh
eval "$(conda shell.bash hook)"
conda activate lpa

# make directories in work
WORKDIR=$WORK/cram_crai
mkdir -p $WORKDIR
cd $WORKDIR

cram_file_path=$1
crai_file_path=$cram_file_path.crai
out_prefix=$(basename $cram_file_path .cram)
out_prefix="${out_prefix}_subset"

python_script_path=$2

conda list | grep pysam

python $python_script_path \
    --input_file $cram_file_path \
    --output_file ${out_prefix}.cram \
    --chrom 6 \
    --threads $SLURM_CPUS_PER_TASK