#!/bin/bash
#SBATCH --job-name=isaac_imputation 
#SBATCH --mail-type=FAIL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32gb
#SBATCH --output=imputation.log
#SBATCH --error=imputation.err

# Stop immediately if any command fails
set -e

# ====================
# CONFIGURATION
# ====================
PROJECT_DIR="/lscratch/GWAS_practice/DATA"
INPUT_DATA="$PROJECT_DIR/cleaned_data/final_qc/wakanda_final_clean"
OUT_DIR="$PROJECT_DIR/cleaned_data/imputation_ready_strict"

mkdir -p "$OUT_DIR"

# =====================
# MODULE LOADING (CRITICAL FOR BGZIP)
# =====================
echo "Cleaning environment..."
module purge || true

echo "Loading modules..."
module load plink/1.9
module load bcftools/1.10.2


# Verify bcftools is actually running
echo "Using bcftools version: $(bcftools --version | head -n 1)"

echo "========================================================"
echo " PREPARING FILES WITH VALIDATION"
echo "========================================================"

# Loop through Chromosomes 1 to 22
for CHR in {1..22}; do
    echo " -> Processing Chromosome $CHR..."
    
    # 1. PLINK: Create uncompressed VCF
    # Using 'vcf-iid' ensures sample IDs are formatted correctly
    plink --bfile "$INPUT_DATA" \
          --chr $CHR \
          --recode vcf-iid \
          --out "$OUT_DIR/wakanda_chr${CHR}" \
          --keep-allele-order \
          --real-ref-alleles \
          --silent

    # 2. COMPRESS: Use bcftools to create BGZF
    # -O z  : Output compressed VCF (BGZF)
    # -o    : Output filename
    # --threads 4 : Use the requested CPUs
    echo "    Compressing..."
    bcftools view "$OUT_DIR/wakanda_chr${CHR}.vcf" \
             -O z \
             --threads 4 \
             -o "$OUT_DIR/wakanda_chr${CHR}.vcf.gz"

    # 3. VERIFY: Try to index the file
    # If the file is NOT valid BGZF, this command will crash the script.
    # This is the guarantee that the file is safe for the server.
    echo "    Verifying (Indexing)..."
    bcftools index -t "$OUT_DIR/wakanda_chr${CHR}.vcf.gz"

    # 4. CLEANUP: Remove the massive text VCF
    rm "$OUT_DIR/wakanda_chr${CHR}.vcf"
    
    echo "    SUCCESS: Chromosome $CHR is valid BGZF."
done

echo "========================================================"
echo " ALL FILES VERIFIED AND READY"
echo " Files are in: $OUT_DIR"
echo "========================================================"








