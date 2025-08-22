#!/bin/bash
#SBATCH --job-name=merge
#SBATCH --output=slurm/merge.out
#SBATCH --error=slurm/merge.err
#SBATCH --time=2-00:00:00
#SBATCH --mem=320G
#SBATCH --partition=general
#SBATCH -n 1
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ztcaterer@colorado.edu

# Load PLINK 2 first for cleaning
module purge
module load plink/2.00a3LM

# Set working directory
WORKDIR=$WORK/merge_files
mkdir -p "$WORKDIR"
cd "$WORKDIR" || exit 1

# Define input prefixes (no extensions)
SOL="$WORK/genotype_harmonizer/harmonized_SOL"
REF="$WORK/vcf2plink/reference_panel"

if [ ! -f harmonized_SOL_cleaned.bed ] && [ ! -f harmonized_SOL_cleaned.bim ] && [ ! -f harmonized_SOL_cleaned.fam ]; then
  # Step 1: Clean SOL dataset
  echo "Step 1: Clean SOL dataset..."
  plink2 --bfile "$SOL" \
        --rm-dup force-first \
        --max-alleles 2 \
        --set-all-var-ids @:#:\$r:\$a \
        --make-bed \
        --out harmonized_SOL_cleaned \
        --threads 4
else 
  echo "✅ SOL dataset already cleaned. Skipping step 1."
fi

if [ ! -f reference_panel_cleaned.bed ] && [ ! -f reference_panel_cleaned.bim ] && [ ! -f reference_panel_cleaned.fam ]; then
  # Step 2: Clean Reference dataset
  echo "Step 2: Clean Reference dataset..."
  plink2 --bfile "$REF" \
       --rm-dup force-first \
       --max-alleles 2 \
       --set-all-var-ids @:#:\$r:\$a \
       --new-id-max-allele-len 20 missing \
       --make-bed \
       --out reference_panel_cleaned \
       --threads 4
else
  echo "✅ Reference dataset already cleaned. Skipping step 2."
fi

# Check if cleaned files exist before merging
if [ ! -f harmonized_SOL_cleaned.bed ] || [ ! -f harmonized_SOL_cleaned.bim ] || [ ! -f harmonized_SOL_cleaned.fam ]; then
  echo "❌ Cleaned SOL files are missing. Exiting."
  exit 1
fi

if [ ! -f reference_panel_cleaned.bed ] || [ ! -f reference_panel_cleaned.bim ] || [ ! -f reference_panel_cleaned.fam ]; then
  echo "❌ Cleaned reference panel files are missing. Exiting."
  exit 1
fi

# Switch to PLINK 1.9 for merging
module load plink/1.90b6.21

# Define base file names
BASEA=harmonized_SOL_cleaned
BASEB=reference_panel_cleaned

# Prepare merge list
echo "$(pwd)/$BASEB" > merge_list.txt

echo "First-pass merge attempt (mode 6 - tolerant)..."
plink --bfile "$BASEA" \
      --merge-list merge_list.txt \
      --merge-mode 5 \
      --make-bed \
      --out tmp_merge \
      --threads 4

# Check if merge was successful or if there are problematic SNPs
if [[ -f tmp_merge-merge.missnp ]]; then
  echo "❌ Merge failed due to strand/allele issues. Found missnp file."
  echo "📝 Removing problematic SNPs and retrying..."
  
  # Remove problematic SNPs from both datasets
  plink --bfile "$BASEA" \
        --exclude tmp_merge-merge.missnp \
        --make-bed \
        --out ${BASEA}_dropbad \
        --threads 4
  
  plink --bfile "$BASEB" \
        --exclude tmp_merge-merge.missnp \
        --make-bed \
        --out ${BASEB}_dropbad \
        --threads 4
  
  # Create new merge list with cleaned datasets
  echo "$(pwd)/${BASEB}_dropbad" > merge_list_final.txt
  
  # Retry merge with cleaned datasets
  echo "🔄 Retrying merge with cleaned datasets..."
  plink --bfile ${BASEA}_dropbad \
        --merge-list merge_list_final.txt \
        --make-bed \
        --out merged_SOL_ref \
        --threads 4

elif [[ -f tmp_merge.bed ]]; then
  echo "✅ First-pass merge successful!"
  # Rename files to final output
  mv tmp_merge.bed merged_SOL_ref.bed
  mv tmp_merge.bim merged_SOL_ref.bim
  mv tmp_merge.fam merged_SOL_ref.fam
  mv tmp_merge.log merged_SOL_ref.log
else
  echo "❌ Unexpected merge failure. Check tmp_merge.log for details."
  exit 1
fi

# Final verification
if [[ -f merged_SOL_ref.bed ]] && [[ -f merged_SOL_ref.bim ]] && [[ -f merged_SOL_ref.fam ]]; then
  echo "🎉 Final merge successful!"
  echo "📊 Output files: merged_SOL_ref.{bed,bim,fam}"
  
  # Print summary statistics
  echo "📈 Summary:"
  echo "   - Merged dataset: $(wc -l < merged_SOL_ref.fam) samples"
  echo "   - Total variants: $(wc -l < merged_SOL_ref.bim) SNPs"
else
  echo "❌ Merge still failed. Check merged_SOL_ref.log for details."
  echo "🔍 Troubleshooting suggestions:"
  echo "   1. Check if input files have compatible formats"
  echo "   2. Verify chromosome coding (1-22, X, Y vs chr1-chr22, chrX, chrY)"
  echo "   3. Check for build/strand mismatches"
  exit 1
fi