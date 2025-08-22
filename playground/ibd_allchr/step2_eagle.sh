#!/bin/bash
#SBATCH --job-name=eagle_chr${CHROM}
#SBATCH --output=slurm/eagle_chr${CHROM}.out
#SBATCH --error=slurm/eagle_chr${CHROM}.err
#SBATCH --time=11-00:00:00
#SBATCH --partition=general
#SBATCH --mem=500g
#SBATCH --cpus-per-task=16
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ztcaterer@colorado.edu

module purge
module load anaconda
module load plink/1.90b6.21
# source /nas/longleaf/apps/anaconda/4.3.0/anaconda/etc/profile.d/conda.sh
conda activate lpa

EAGLE="$SOFTWARE/Eagle_v2.4.1/eagle"
MAP="$SOFTWARE/Eagle_v2.4.1/tables/genetic_map_1cMperMb.txt"

cd "$WORK/ibd_allchr/chr${CHROM}"

plink \
 --bfile chr${CHROM} \
 --chr $CHROM \
 --set-missing-var-ids @:# \
 --make-bed \
 --out chr${CHROM}_filtered

"$EAGLE" \
 --numThreads 16 \
 --bfile=chr${CHROM}_filtered \
 --chrom=$CHROM \
 --geneticMapFile="$MAP" \
 --outPrefix=phased_chr${CHROM}
