#!/bin/bash
#SBATCH --job-name=pca
#SBATCH --output=slurm/pca.out
#SBATCH --error=slurm/pca.err
#SBATCH --time=10-00:00:00
#SBATCH --mem=100G
#SBATCH --partition=general
#SBATCH -n 1
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ztcaterer@colorado.edu

# Load PLINK 2
module purge
module load plink/2.00

# Set working directory and input
WORKPCA=$WORK/pca
mkdir -p "$WORKPCA"
cd $WORKPCA || exit 1
MERGE_SOL_INPUT="$WORK/merge_files/merged_SOL_ref"
OUT=merged_SOL_ref

INPUT="merged_SOL_ref"
if [ ! -f "${INPUT}.bed" ] || [ ! -f "${INPUT}.bim" ] || [ ! -f "${INPUT}.fam" ]; then
    echo "Input files not found in $PWD. Copying from $MERGE_SOL_INPUT..."
    cp "${MERGE_SOL_INPUT}.bed" .
    cp "${MERGE_SOL_INPUT}.bim" .
    cp "${MERGE_SOL_INPUT}.fam" .
fi

if [ ! -f "${OUT}_maf05.bed" ] || [ ! -f "${OUT}_maf05.bim" ] || [ ! -f "${OUT}_maf05.fam" ]; then
    # Run MAF filtering with threshold 0.05
    echo "Running MAF filtering..."
    plink2 --bfile $INPUT \
      --maf 0.05 \
      --make-bed \
      --out ${OUT}_maf05
fi

if [ ! -f "${OUT}_pruned.prune.in" ]; then
    echo "LD pruning not found. Running LD pruning..."
    
    # Run LD pruning (window size = 50 variants, step size = 5, r^2 threshold = 0.2)
    echo "Running LD pruning..."
    plink2 --bfile ${OUT}_maf05 \
      --indep-pairwise 50 5 0.2 \
      --out ${OUT}_pruned

    echo "LD pruning complete. Output saved to ${OUT}_pruned.prune.in"
    echo "Extracting pruned SNPs..."
    plink2 --bfile ${OUT}_maf05 \
        --extract ${OUT}_pruned.prune.in \
        --make-bed \
        --out ${OUT}_pruned
        
else
    echo "LD pruning already done."
fi

# Ensure the output files are in the working directory
if [ -f "${OUT}.eigenvec" ] && [ -f "${OUT}.eigenval" ]; then
    echo "PCA output files are ready in $WORKPCA"
else
    echo "Running PCA..."

    # Run PCA on LD-pruned and MAF-filtered data
    plink2 --bfile ${OUT}_maf05 \
    --extract ${OUT}_pruned.prune.in \
    --pca 20 \
    --threads 4 \
    --out $OUT

    echo "✅ PCA complete. Output saved to $OUT.eigenvec and $OUT.eigenval"
fi

# Run R script for plotting
module load r

if [ ! -f "plot_pca.R" ]; then
    cp $LPA/playground/local_ancestry/plot_pca.R .
fi

if [ ! -f "figures/" ]; then
    mkdir figures
fi

for i in {1..9}; do
    Rscript plot_pca.R $i $((i+1))
done

Rscript plot_pca.R 1 3


