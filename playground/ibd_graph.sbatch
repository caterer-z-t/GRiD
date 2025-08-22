#!/bin/bash
#SBATCH --job-name=ibd_graph_spring
#SBATCH --output=slurm/ibd_graph_spring.out
#SBATCH --error=slurm/ibd_graph_spring.err
#SBATCH --time=10:00:00
#SBATCH --partition=general
#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH --mem=100g
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ztcaterer@colorado.edu

# Exit on error and undefined vars
set -euo pipefail

module purge
module load anaconda/4.3.0
source /nas/longleaf/apps/anaconda/4.3.0/anaconda/etc/profile.d/conda.sh
conda activate lpa

timestamp() {
  date +"[%Y-%m-%d %H:%M:%S]"
}
log() {
  echo "$(timestamp) $1"
}
log "🚀 Starting IBD Figure Script"

conda update scipy
# or
# pip install --upgrade scipy

python "$LPA/playground/ibd_graph.py"
