#!/bin/bash
#SBATCH --job-name=compare_plink
#SBATCH --output=slurm/compare_plink.out
#SBATCH --error=slurm/compare_plink.err
#SBATCH --time=02:00:00
#SBATCH --partition=general
#SBATCH -n 1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32g
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ztcaterer@colorado.edu

# Exit on errors
set -euo pipefail
echo "Job started on $(hostname) at $(date)"

# Load PLINK 1.90 (supports --bmerge)
module purge
module load plink/1.90b6.21

# Define input prefixes
SET_A_PREFIX=$LPA/vcf_to_plink/plink_output/cardiac_chr6
SET_B_PREFIX=$LPA/vcf_to_plink/plink_output_cleaned/cardiac_chr6

WORKDIR=$WORK/compare_plink_chr6
mkdir -p "$WORKDIR"
cd "$WORKDIR"

# Compare .fam sample IDs
echo "Comparing .fam sample IDs..."
cut -d' ' -f1-2 "${SET_A_PREFIX}.fam" > samples_A.txt
cut -d' ' -f1-2 "${SET_B_PREFIX}.fam" > samples_B.txt
diff samples_A.txt samples_B.txt > fam_diff.txt || echo "✅ Sample ID differences found (saved in fam_diff.txt)"

# Compare .bim variant IDs
echo "Comparing .bim variant IDs..."
cut -f2 "${SET_A_PREFIX}.bim" > snps_A.txt
cut -f2 "${SET_B_PREFIX}.bim" > snps_B.txt
comm -3 <(sort snps_A.txt) <(sort snps_B.txt) > bim_diff.txt || echo "✅ Variant ID differences found (saved in bim_diff.txt)"

# Try to merge using --merge-mode 7 to compare genotypes
echo "Comparing genotypes using PLINK --merge-mode 7..."
plink \
  --bfile "$SET_A_PREFIX" \
  --bmerge "${SET_B_PREFIX}.bed" "${SET_B_PREFIX}.bim" "${SET_B_PREFIX}.fam" \
  --merge-mode 7 \
  --out merge_compare || echo "⚠️ Merge warnings encountered. See merge_compare.log"

# Report results
if [ -s merge_compare.diff ]; then
  echo "❌ Genotype differences found! See merge_compare.diff"
else
  echo "✅ No genotype differences found."
fi

echo "✅ PLINK file comparison completed at $(date)"

# Clean up temporary files
cd ..
rm -rf "$WORKDIR"
echo "Temporary files cleaned up."
echo "Job completed at $(date)"
