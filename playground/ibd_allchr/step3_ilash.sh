#!/bin/bash
#SBATCH --job-name=ilash_chr${CHROM}
#SBATCH --output=slurm/ilash_chr${CHROM}.out
#SBATCH --error=slurm/ilash_chr${CHROM}.err
#SBATCH --time=11-00:00:00
#SBATCH --partition=general
#SBATCH --mem=500g
#SBATCH --cpus-per-task=32
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ztcaterer@colorado.edu

module purge
module load ilash/1.0.2
module load anaconda
# source /nas/longleaf/apps/anaconda/4.3.0/anaconda/etc/profile.d/conda.sh
conda activate lpa

cd "$WORK/ibd_allchr/chr${CHROM}"

gunzip -c phased_chr${CHROM}.haps.gz > phased_chr${CHROM}.haps

python \
 "$SOFTWARE/ilash_analyzer/zc_hap_sample_to_ped.py" \
 "phased_chr${CHROM}" \
 "$SOFTWARE/Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz" \
 "$CHROM" \
 "$(wc -l < phased_chr${CHROM}.haps)" \
 "phased_chr${CHROM}"

cat > config_chr${CHROM}.txt <<EOF
map phased_chr${CHROM}.map
ped phased_chr${CHROM}.ped
output output_chr${CHROM}.match
slice_size 350
step_size 350
perm_count 20
shingle_size 15
shingle_overlap 0
bucket_count 5
max_thread 32
match_threshold 0.99
interest_threshold 0.70
min_length 2.9
auto_slice 1
slice_length 2.9
cm_overlap 1
minhash_threshold 55
EOF

ilash config_chr${CHROM}.txt
