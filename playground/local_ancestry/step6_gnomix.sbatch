#!/bin/bash
#SBATCH --job-name=gnomix
#SBATCH --output=slurm/%x.out
#SBATCH --error=slurm/%x.err
#SBATCH --time=11-00:00:00
#SBATCH --mem=300G
#SBATCH -n 1
#SBATCH --cpus-per-task=16
#SBATCH --mail-type=ALL 
#SBATCH --mail-user=ztcaterer@colorado.edu

set -euo pipefail

###############################
# Load modules
###############################

module purge
module load anaconda 
source $(conda info --base)/etc/profile.d/conda.sh
conda activate gnomix
module load bcftools
module load samtools
module load plink

###############################
# Paths and variables
###############################
export NUMEXPR_MAX_THREADS=$SLURM_CPUS_PER_TASK

WORKDIR=$WORK/gnomix_local_ancestry
FILEDIR=$WORKDIR/files
OUTDIR=$WORKDIR/output
mkdir -p $WORKDIR $FILEDIR $OUTDIR

cd $WORKDIR

GNOMIX_PY=$WORKDIR/gnomix/gnomix.py

IN_VCF_BASE="${WORK}/phasing/phased_merged_SOL_ref"
IN_VCF="${IN_VCF_BASE}.vcf"
IN_VCF_GZ="${IN_VCF_BASE}.vcf.gz"

# Ensure compress vcf file to gunzip
bgzip -c $IN_VCF > $IN_VCF_GZ
tabix -p vcf $IN_VCF_GZ

###############################
# Copy gnomix to work dir
###############################
if [ ! -f "$GNOMIX_PY" ]; then
    cp -rf $SOFTWARE/gnomix $WORKDIR/
fi

exit 1 # exit early for now to debug


###############################
# Split Phased VCF into reference panel + sample subsets
###############################

# for any subject id that matches to subject id in ./files/gnomad_meta_key.txt, pull out 
# of phased VCF and put into phased reference VCF 
# Prepare clean list of reference samples for bcftools
REF_IDS="${GRID}/playground/local_ancestry/files/gnomad_meta_key.txt"
REF_IDS_NOPAD="${FILEDIR}/gnomad_meta_key_ids.txt"
awk 'NR>1 {print $2}' $REF_IDS > $REF_IDS_NOPAD  # assuming IID matches VCF sample names

bcftools view -S $REF_IDS_NOPAD --force-samples -Oz -o ${FILEDIR}/phased_ref.vcf.gz $IN_VCF_GZ
bcftools index ${FILEDIR}/phased_ref.vcf.gz


# now lets make the genetic map from VCF file using the python file provided by jonathan
MAKE_GEN_MAP=$GRID/playground/local_ancestry/files/make_genetic_map.py
python $MAKE_GEN_MAP \
    --vcf ${FILEDIR}/phased_ref.vcf.gz \
    --map $SOFTWARE/Eagle_v2.4.1/tables/genetic_map_hg38_chr6.txt.gz \
    --out ${FILEDIR}/genetic_map.tsv \
    --chrom 6

# now lets make the sample map file matching reference samples to their respective reference populations
SAMPLE_MAP=$GRID/playground/local_ancestry/files/gnomad_sample_map.txt
cp $SAMPLE_MAP ${FILEDIR}/sample_map.txt

# filter this sample map to only include samples present in the phased_ref.vcf.gz
bcftools query -l ${FILEDIR}/phased_ref.vcf.gz > ${FILEDIR}/ref_sample_ids.txt
grep -Ff ${FILEDIR}/ref_sample_ids.txt ${FILEDIR}/sample_map.txt > ${FILEDIR}/sample_map_filtered.txt
mv ${FILEDIR}/sample_map_filtered.txt ${FILEDIR}/sample_map.txt

VCF_OUT=${FILEDIR}/phased_samples_gnomix.vcf.gz

echo "🧬 Running gnomix on ${VCF_OUT}"
python $GNOMIX \
    "$VCF_OUT" \
    "$OUTDIR" \
    6 \
    True \
    ${FILEDIR}/genetic_map.tsv \
    ${FILEDIR}/phased_ref.vcf.gz \
    ${FILEDIR}/sample_map.txt \
    
echo "✅ gnomix completed at $(date)"