#!/bin/bash
#SBATCH --job-name=convert_to_ped
#SBATCH --output=slurm/convert_to_ped.out
#SBATCH --error=slurm/convert_to_ped.err
#SBATCH --time=10-00:00:00
#SBATCH --partition=general
#SBATCH -n 1
#SBATCH --cpus-per-task=10
#SBATCH --mem=900g
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ztcaterer@colorado.edu

set -euo pipefail
echo "Job started on $(hostname) at $(date)"

# Load modules
module purge
# module load plink/1.90b6.21
# module load shapeit/2.837
module load anaconda/4.3.0
source /nas/longleaf/apps/anaconda/4.3.0/anaconda/etc/profile.d/conda.sh
conda activate lpa

# Define paths (edit these as needed)
LPA=/nas/longleaf/home/catererz/epi/lpa
WORKDIR=$WORK/convert_to_ped

mkdir -p "$WORKDIR"
cd "$WORKDIR"

HAPSFILE=$LPA/eagle_output_files/lpa_gene/phased_chr6_LPA_plink.haps
HAPSFILEGUNZIP=$HAPSFILE.gz
SAMPLEFILE=$LPA/eagle_output_files/lpa_gene/phased_chr6_LPA_plink.sample

cp "$HAPSFILE" .
cp "$SAMPLEFILE" .
cp "$HAPSFILEGUNZIP" .

HAPS_BASENAME=$(basename "$HAPSFILE" .haps)

haps_sample_to_ped_script="$LPA/playground/zc_hap_sample_to_ped.py"

genetic_map_file="$SOFTWARE/Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz"
# Get SNP count from haps file
snp_count=$(wc -l < "${HAPS_BASENAME}.haps")
echo "Number of SNPs in haps file: $snp_count"

python "$haps_sample_to_ped_script" \
    $HAPS_BASENAME \
    $genetic_map_file \
    6 \
    $snp_count \
    output_phased \
    
# Sanity checks
if [[ -f output_phased.ped ]]; then
    echo "✅ output_phased.ped created successfully."
else
    echo "❌ output_phased.ped does not exist!" >&2
    exit 1
fi

# Peek at contents
head -2 output_phased.ped
head -5 "${HAPS_BASENAME}.haps"

# Move output to final destination
OUTPUT_DIR=$LPA/playground/convert_to_ped
mkdir -p "$OUTPUT_DIR"

mv output_phased.ped "$OUTPUT_DIR/lpa_gene_phased.ped"
mv output_phased.map "$OUTPUT_DIR/lpa_gene_phased.map"

# Clean up
cd ../
rm -rf "$WORKDIR"

echo "✅ Job completed at $(date)"
