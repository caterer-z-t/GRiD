#!/bin/bash
#SBATCH --job-name=vcf2plink_chr${CHROM}
#SBATCH --output=slurm/vcf2plink_chr${CHROM}.out
#SBATCH --error=slurm/vcf2plink_chr${CHROM}.err
#SBATCH --time=11-00:00:00
#SBATCH --partition=general
#SBATCH --mem=500g
#SBATCH --cpus-per-task=12
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ztcaterer@colorado.edu

module purge
module load plink/2.00a
module load samtools/1.21

VCF="$SOL/sequenced/cardiac_cohorts_SOL_dragen.chr${CHROM}.vcf.gz"
WORKDIR="$WORK/ibd_allchr/chr${CHROM}"
mkdir -p "$WORKDIR"
cd "$WORKDIR"

plink2 \
 --vcf "$VCF" \
 --make-bed \
 --out chr${CHROM} \
 --threads 12
