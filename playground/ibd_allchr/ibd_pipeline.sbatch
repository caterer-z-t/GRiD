#!/bin/bash
#SBATCH --job-name=ibd_pipeline
#SBATCH --output=slurm/ibd_pipeline.out
#SBATCH --error=slurm/ibd_pipeline.err
#SBATCH --time=11-00:00:00
#SBATCH --partition=general
#SBATCH -n 1
#SBATCH --cpus-per-task=32
#SBATCH --mem=500g
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ztcaterer@colorado.edu

# Exit on error and undefined vars
set -euo pipefail

timestamp() {
  date +"[%Y-%m-%d %H:%M:%S]"
}

log() {
  echo "$(timestamp) $1"
}

log "🚀 Starting IBD pipeline"

# Load modules
module purge
module load plink/2.00a
module load samtools/1.21

# Define paths
WORKDIR="$WORK/ibd_pipeline"
VCF_SOURCE="$SOL/sequenced/cardiac_cohorts_SOL_dragen.chr6.vcf.gz"
VCF_NAME="$(basename "$VCF_SOURCE")"
CLEANED_VCF="cleaned_${VCF_NAME}"
OUTPUT_DIR="$LPA/ibd_pipeline"

mkdir -p "$WORKDIR" "$OUTPUT_DIR"
cd "$WORKDIR"
log "📁 Working in $WORKDIR"

# Copy and clean VCF
if [[ ! -f "$VCF_NAME" ]]; then
  log "📄 Copying VCF file"
  cp "$VCF_SOURCE" .
else
  log "📄 VCF file already exists, skipping copy"
fi

if [[ ! -f "$CLEANED_VCF" ]]; then
  log "🧹 Cleaning VCF (removing malformed line 2336741)"
  zcat "$VCF_NAME" | sed '2336741d' | bgzip -c > "$CLEANED_VCF"
  tabix -p vcf "$CLEANED_VCF"
else
  log "✅ Cleaned VCF already exists, skipping cleanup"
fi

# Run PLINK
if [[ ! -f "chr6.bed" ]]; then
  log "🔧 Running PLINK on cleaned VCF"
  plink2 --vcf "$CLEANED_VCF" --make-bed --out chr6 --threads 32
else
  log "✅ PLINK outputs exist, skipping"
fi

# Sanity check
for f in chr6.{bed,bim,fam}; do
  [[ -f $f ]] || { log "❌ Missing $f"; exit 1; }
done

# Clean intermediate files
rm -f "$CLEANED_VCF" "$VCF_NAME"
log "🧼 Cleaned up VCF intermediate files"

# Step 2: Eagle phasing
module purge
module load anaconda/4.3.0
module load plink/1.90b6.21
source /nas/longleaf/apps/anaconda/4.3.0/anaconda/etc/profile.d/conda.sh
conda activate lpa

EAGLE_EXE="$SOFTWARE/Eagle_v2.4.1/eagle"
GEN_MAP="$SOFTWARE/Eagle_v2.4.1/tables/genetic_map_1cMperMb.txt"

cp "$GEN_MAP" .

log "📍 Filtering chr6 LPA region and assigning SNP IDs"
plink --bfile chr6 --chr 6 --from-bp 110000000 --to-bp 210000000 \
  --set-missing-var-ids @:# --make-bed --out chr6_LPA

log "🧬 Running Eagle"
"$EAGLE_EXE" --numThreads 32 --bfile=chr6_LPA --chrom=6 \
  --geneticMapFile="$(basename "$GEN_MAP")" \
  --outPrefix=phased_chr6_LPA

[[ -f phased_chr6_LPA.haps.gz && -f phased_chr6_LPA.sample ]] || {
  log "❌ Eagle output missing!"; exit 1;
}

rm -f chr6_LPA.{bed,bim,fam}
log "✅ Eagle LPA phasing completed"

# Step 3: Convert phased output to PLINK PED/MAP
module purge
module load anaconda/4.3.0
source /nas/longleaf/apps/anaconda/4.3.0/anaconda/etc/profile.d/conda.sh
conda activate lpa

HAPSFILEGUNZIP="phased_chr6_LPA.haps.gz"
HAPS_BASENAME="phased_chr6_LPA"
haps_sample_to_ped_script="$SOFTWARE/ilash_analyzer/zc_hap_sample_to_ped.py"
genetic_map_file="$SOFTWARE/Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz"

log "📤 Unzipping haps file for processing"
gunzip -c "$HAPSFILEGUNZIP" > "${HAPS_BASENAME}.haps"
snp_count=$(wc -l < "${HAPS_BASENAME}.haps")
log "🔢 SNP count = $snp_count"

python "$haps_sample_to_ped_script" \
  "$HAPS_BASENAME" "$genetic_map_file" 6 "$snp_count" "$HAPS_BASENAME"

for ext in ped map; do
  [[ -f "${HAPS_BASENAME}.${ext}" ]] || { log "❌ Missing .${ext} file!"; exit 1; }
done

log "✅ Step 3 completed. Generated .ped and .map files"

# Step 4: Run iLASH
module purge
module load gcc/12.2.0
module load ilash/1.0.2
module load anaconda/4.3.0
source /nas/longleaf/apps/anaconda/4.3.0/anaconda/etc/profile.d/conda.sh
conda activate lpa

cat > config.txt <<EOF
map ${HAPS_BASENAME}.map
ped ${HAPS_BASENAME}.ped
output output.match
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

log "🚀 Running iLASH"
ilash config.txt

[[ -f output.match ]] || { log "❌ iLASH output.match not found!"; exit 1; }
log "✅ iLASH completed"

# Move and cleanup
log "📦 Moving results to $OUTPUT_DIR"
mv -f ./* "$OUTPUT_DIR"

cd ../ && rm -rf "$WORKDIR"
log "✅ All steps completed successfully"
