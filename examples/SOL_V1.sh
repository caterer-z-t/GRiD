#!/bin/bash
#SBATCH --job-name=SOL_PIPELINE
#SBATCH --output=slurm/%x.out
#SBATCH --error=slurm/%x.err
#SBATCH --time=2-00:00:00
#SBATCH --partition=general
#SBATCH -n 1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ztcaterer@colorado.edu

# ------------------------------
# Environment Setup
# ------------------------------
module purge
module load anaconda
source $(conda info --base)/etc/profile.d/conda.sh 
eval "$(conda shell.bash hook)" 
conda activate GRiD

export WORK=/work/users/c/a/catererz
export USERS=/users/c/a/catererz
export HOME=/nas/longleaf/home/catererz
export EPI=/proj/epi
export SOL=$EPI/Genetic_Data_Center/SOL
export CRAM=$SOL/cram
export LPA=$HOME/epi/lpa
export GRID=$HOME/epi/GRiD
export SOFTWARE=$HOME/epi/software

# Name for working directory for this run
WORK_DIR=$WORK/software_dev/sol_pipeline_run

mkdir -p $WORK_DIR
cd $WORK_DIR

grid wgs $GRID/grid/example_config.yaml