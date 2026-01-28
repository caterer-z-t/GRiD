#!/bin/bash
#SBATCH --job-name=find_freeze
#SBATCH --output=slurm/%x.out
#SBATCH --error=slurm/%x.err
#SBATCH --time=10-00:00:00
#SBATCH --partition=general
#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ztcaterer@colorado.edu

module purge

find / -name "hchs_genomics_id_mapping_240715.txt" 