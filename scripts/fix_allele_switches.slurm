#!/bin/bash
#SBATCH --job-name=Fix_Alleles_switches
#SBATCH --mail-type=FAIL
#SBATCH --output=fix_alleles_switches.log
#SBATCH --mem=16gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4

set -e

# ==============================================================================
# CONFIGURATION
# ==============================================================================
# Input: The files that were verified
INPUT_DIR="cleaned_data/imputation_ready_strict"

# Output: The new, corrected files
OUT_DIR="cleaned_data/imputation_ready_strict_fix"
mkdir -p "$OUT_DIR"

# Reference Genome (The one you just downloaded)
REF_GENOME="/home/abegunde/DATA_OMICS_PROJECT/DATA/human_g1k_v37.fasta"

# ==============================================================================
# MODULES
# ==============================================================================
module purge || true
module load bcftools/1.10.2

echo "Processing Chromosomes..."

for CHR in {1..22}; do
    echo " -> Correcting Chromosome $CHR..."
    
    # The Magic Command: bcftools norm
    # --check-ref s : "Set" the Ref allele to match the genome (Swaps alleles if needed)
    # -f : The reference genome file
    
    bcftools norm \
        -f "$REF_GENOME" \
        --check-ref s \
        --threads 4 \
        -O z \
        -o "$OUT_DIR/wakanda_chr${CHR}.vcf.gz" \
        "$INPUT_DIR/wakanda_chr${CHR}.vcf.gz"

    # Index the new file
    bcftools index -t "$OUT_DIR/wakanda_chr${CHR}.vcf.gz"
    
    echo "    Fixed Chromosome $CHR."
done

echo "Done! Fixed files are in $OUT_DIR"
