#!/bin/bash
#SBATCH --job-name=ilash
#SBATCH --output=slurm/ilash.out
#SBATCH --error=slurm/ilash.err
#SBATCH --time=11-00:00:00
#SBATCH --mem=300G
#SBATCH -n 1
#SBATCH --cpus-per-task=32
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ztcaterer@colorado.edu


set -euo pipefail
echo "Job started on $(hostname) at $(date)"

# Load modules
module purge
module load anaconda
module load bcftools
module load samtools

# Initialize Conda
source $(conda info --base)/etc/profile.d/conda.sh
conda activate gnomix

# Variables
WORKDIR="$WORK/ibd"
FILEDIR="$WORK/gnomix/output/"
SOFTWAREDIR="$SOFTWARE/iLASH"
ILASH="$SOFTWAREDIR/ilash"
CHROM=6  # or dynamically pass this via sbatch --export=CHROM=6

# Ensure working directory exists
mkdir -p "$WORKDIR"
cd "$WORKDIR"

# Copy split files if they don't exist
if ls gnomix*[0-9].vcf 1> /dev/null 2>&1; then
    echo "Split VCF files already exist, skipping copy."
else
    echo "Copying split VCF files..."
    cp "$FILEDIR"/*.vcf "$WORKDIR/"
fi

echo "Compressing and indexing VCF files..."

for vcf in *.vcf; do
    if [[ ! -f "${vcf}.gz" ]]; then
        echo "Compressing $vcf..."
        bgzip "$vcf"
        tabix -p vcf "${vcf}.gz"
    else
        echo "$vcf already compressed, skipping."
    fi
done

# Merge all VCFs into one big VCF first
MERGED_VCF="phased_chr${CHROM}.vcf"
if [ -f "$MERGED_VCF" ]; then
    echo "Merged VCF already exists, skipping."
else
    echo "Merging VCFs into $MERGED_VCF..."
    bcftools merge -m none -Ov -o "$MERGED_VCF" gnomix*.vcf.gz
fi

# Count SNPs
SNP_COUNT=$(bcftools view "$MERGED_VCF" | grep -v '^#' | wc -l)
echo "SNP count in merged VCF: $SNP_COUNT"

# Convert merged VCF to PED
VCFTOPED="$SOFTWARE/ilash_analyzer/vcf2ped.py"
BASENAME=$(basename "$MERGED_VCF" .vcf)
PED_OUT="${BASENAME}"

# vcf2ped.py command line arguments
# python <vcf-address> \
#     <output_address> \
#     <snp_count>

echo "Running vcf2ped..."
# python "$VCFTOPED" "$MERGED_VCF" "$PED_OUT" "$SNP_COUNT"

# Generate map file
MAP_FILE="phased_chr${CHROM}.map"
if [ -f "$MAP_FILE" ]; then
    echo "Map file already exists."
else
    echo "Generating map file..."
    bcftools query -f '%CHROM\t%ID\t0\t%POS\n' "$MERGED_VCF" > "$MAP_FILE"
    bcftools query -f '%CHROM\t%POS\t%ID[\t%GT]\n' "$MERGED_VCF"

fi


# iLASH config file
CONFIG="ilash_config_chr${CHROM}.txt"
if [ -f "$CONFIG" ]; then
    echo "iLASH config file already exists, skipping creation."
else
    echo "Creating iLASH config file..."
    cat > "$CONFIG" <<EOF
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
fi

# Run iLASH
if [ -x "$ILASH" ]; then
    echo "Running iLASH..."
    "$ILASH" "$CONFIG"
else
    echo "ERROR: iLASH binary not found or not executable at $ILASH"
    exit 1
fi

echo "Job finished at $(date)"
