#!/bin/bash
#SBATCH --job-name=GWAS_combined
#SBATCH --mem=16gb
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/GWAS_combined_%j.log

# ================= CONFIGURATION =================
INPUT_DIR="results/final_gwas_results"
OUTPUT_DIR="results/summary_statistics"

mkdir -p "$OUTPUT_DIR"

echo "Starting Result Collection..."

# 1. Detect Phenotypes automatically
# We look at Chromosome 1 files to see what traits were analyzed.
# Filename format: chr1_results.[PHENOTYPE].assoc.linear
PHENO_LIST=$(ls $INPUT_DIR/chr1_results.*.assoc.linear | sed 's/.*chr1_results\.//' | sed 's/\.assoc\.linear//')

if [ -z "$PHENO_LIST" ]; then
    echo "CRITICAL ERROR: No result files found in $INPUT_DIR!"
    echo "Did the GWAS jobs finish?"
    exit 1
fi

echo "Detected Phenotypes: "
echo "$PHENO_LIST"
echo "========================================================"

# 2. Loop through each phenotype and merge
for PHENO in $PHENO_LIST; do
    echo "Processing: $PHENO"
    
    FINAL_FILE="${OUTPUT_DIR}/${PHENO}_GWAS_SUMMARY.txt"
    TOP_HITS="${OUTPUT_DIR}/${PHENO}_TOP_HITS.txt"
    
    # --- Step A: Initialize the file with a Header ---
    # We grab the header from Chr 1
    head -n 1 "$INPUT_DIR/chr1_results.${PHENO}.assoc.linear" > "$FINAL_FILE"
    
    # --- Step B: Concatenate all 22 chromosomes ---
    # We assume file names are consistent. We skip the header (tail -n +2)
    for i in {1..22}; do
        FILE="$INPUT_DIR/chr${i}_results.${PHENO}.assoc.linear"
        
        if [ -f "$FILE" ]; then
            tail -n +2 "$FILE" >> "$FINAL_FILE"
        else
            echo "WARNING: Missing results for Chr $i (Pheno: $PHENO)"
        fi
    done
    
    # --- Step C: Extract Significant Hits (P < 5e-8) ---
    # awk column 9 is usually 'P' in PLINK linear output, but let's be safe.
    # We filter for scientific notation (e.g. 1.2e-09) or standard low numbers.
    
    # Get Header
    head -n 1 "$FINAL_FILE" > "$TOP_HITS"
    # Filter rows where P-value (Col 9) is less than 5e-8
    awk '$9 < 0.00000005' "$FINAL_FILE" >> "$TOP_HITS"
    
    HIT_COUNT=$(wc -l < "$TOP_HITS")
    ((HIT_COUNT--)) # Subtract header
    
    echo "  -> Saved Full Summary: $FINAL_FILE"
    if [ "$HIT_COUNT" -gt 0 ]; then
        echo "  -> *** FOUND $HIT_COUNT SIGNIFICANT HITS! *** (See $TOP_HITS)"
    else
        echo "  -> No genome-wide significant hits found."
    fi
    echo "--------------------------------------------------------"

done

echo "Collection Complete!"
