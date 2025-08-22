#!/bin/bash
#SBATCH --job-name=unzip
#SBATCH --output=slurm/unzip.out
#SBATCH --error=slurm/unzip.err
#SBATCH --time=10:00:00
#SBATCH --partition=general
#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH --mem=320g
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ztcaterer@colorado.edu

set -euo pipefail
echo "Job started on $(hostname) at $(date)"

module purge

HAPS_GZ=/nas/longleaf/home/catererz/epi/lpa/eagle_output_files/ch6/phased_chr6_plink.haps.gz
HAPS_OUT=/nas/longleaf/home/catererz/epi/lpa/eagle_output_files/ch6/phased_chr6_plink.haps

# Remove if exists to prevent stale handle issues
if [[ -f "$HAPS_OUT" ]]; then
    echo "Removing stale file $HAPS_OUT"
    rm "$HAPS_OUT"
fi

# Unzip safely
gunzip -c "$HAPS_GZ" > "$HAPS_OUT"

echo "✅ Unzipping complete at $(date)"
