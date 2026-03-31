#!/bin/bash
#SBATCH --job-name=isaac_Merge
#SBATCH --mem=64gb
#SBATCH --cpus-per-task=4
#SBATCH --output=logs/Merge_chromosome.log

# ================= CONFIGURATION =================
# INPUT 
CLEAN_DIR="cleaned_data/imputation_ready_strict_fix/michigan_result/imputation_clean"
# OUTPUT
OUT_DIR="cleaned_data/imputation_ready_strict_fix/final_merged_data"

module load plink/1.9 || module load plink

mkdir -p "$OUT_DIR"

# ================= STEP 1: MAKE THE LIST =================
echo "Generating list of files to merge..."
LIST_FILE="${OUT_DIR}/merge_list.txt"
rm -f "$LIST_FILE"

for i in {1..22}; do
    FILE="${CLEAN_DIR}/chr${i}_imputed_clean"
    if [ -f "${FILE}.bed" ]; then
        echo "$FILE" >> "$LIST_FILE"
    else
        echo "ERROR: Missing Chromosome $i"
    fi
done

# ================= STEP 2: ATTEMPT MERGE =================
echo "Attempting Merge..."

plink \
    --merge-list "$LIST_FILE" \
    --make-bed \
    --allow-no-sex \
    --memory 60000 \
    --out "${OUT_DIR}/wakanda_imputed_final"

# ================= STEP 3: THE RESCUE (If Step 2 Fails) =================
# If PLINK fails due to conflicts, it generates a 'missnp' file.
if [ -f "${OUT_DIR}/wakanda_imputed_final-merge.missnp" ]; then
    
    echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    echo "MERGE CONFLICT DETECTED! (Multiallelic Variants)"
    echo "Starting Automatic Rescue..."
    echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"

    # 1. Flip the bad variants in every chromosome
        
    BAD_SNPS="${OUT_DIR}/wakanda_imputed_final-merge.missnp"
    mkdir -p "${OUT_DIR}/temp_fix"
    NEW_LIST_FILE="${OUT_DIR}/merge_list_fixed.txt"
    rm -f "$NEW_LIST_FILE"

    # Loop through all 22 chromosomes again to remove the bad SNPs
    for i in {1..22}; do
        echo "Excluding bad SNPs from Chr $i..."
        
        plink \
            --bfile "${CLEAN_DIR}/chr${i}_imputed_clean" \
            --exclude "$BAD_SNPS" \
            --make-bed \
            --out "${OUT_DIR}/temp_fix/chr${i}_fixed"
            
        echo "${OUT_DIR}/temp_fix/chr${i}_fixed" >> "$NEW_LIST_FILE"
    done

    # 2. Try Merging Again with the fixed files
    echo "Re-attempting Merge with fixed files..."
    
    plink \
        --merge-list "$NEW_LIST_FILE" \
        --make-bed \
        --allow-no-sex \
        --memory 60000 \
        --out "${OUT_DIR}/wakanda_imputed_final"

    echo "Rescue Complete. Check if 'wakanda_imputed_final.bed' exists."
    
     cleanup of temp files 
    rm -rf "${OUT_DIR}/temp_fix"

else
    echo "Merge was successful on the first try!"
fi
