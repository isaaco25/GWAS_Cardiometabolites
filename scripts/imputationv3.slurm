#!/bin/bash
#SBATCH --job-name=Wakanda_Imp_Prep
#SBATCH --mail-type=FAIL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8gb
#SBATCH --time=04:00:00
#SBATCH --output=imputation_%j.log
#SBATCH --error=imputation_%j.err

# STOP SCRIPT ON ERROR (Prevents false "Success" messages)
set -e

# ==============================================================================
# CONFIGURATION
# ==============================================================================
PROJECT_DIR="/home/abegunde/DATA_OMICS_PROJECT/DATA"
INPUT_DATA="$PROJECT_DIR/cleaned_data/final_qc/wakanda_final_clean"
OUT_DIR="$PROJECT_DIR/cleaned_data/imputation_ready"

mkdir -p "$OUT_DIR"

# ==============================================================================
# MODULE SETUP
# ==============================================================================
echo "Cleaning environment..."
# Ignore errors during purge 
module purge || true

echo "Loading specific modules..."
module load plink/1.9
module load bcftools/1.10.2
module load htslib/1.9

# ==============================================================================
# EXECUTION
# ==============================================================================

# Loop through Chromosomes 1 to 22
for CHR in {1..22}; do
    echo " -> Processing Chromosome $CHR..."
    
    # 1. PLINK: Convert to VCF
    # We proceed because bcftools usually handles the resulting VCF fine.
    plink --bfile "$INPUT_DATA" \
          --chr $CHR \
          --recode vcf-iid \
          --out "$OUT_DIR/wakanda_chr${CHR}" \
          --keep-allele-order \
          --real-ref-alleles \
          --silent

    # 2. BCFTOOLS: Compress to .vcf.gz (BGZF)
    # FIX: Used -O z (Output Type) and -o (Output File) to avoid ambiguity
    bcftools view "$OUT_DIR/wakanda_chr${CHR}.vcf" \
             -O z \
             --threads 4 \
             -o "$OUT_DIR/wakanda_chr${CHR}.vcf.gz"

    # 3. CLEANUP
    rm "$OUT_DIR/wakanda_chr${CHR}.vcf"
    
    echo "    Created: $OUT_DIR/wakanda_chr${CHR}.vcf.gz"
done

echo "========================================================"
echo " PREPARATION COMPLETE"
echo "========================================================"
