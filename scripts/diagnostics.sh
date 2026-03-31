#!/bin/bash
# run_diagnostics.sh

# CONFIG
PROJECT_DIR="/home/abegunde/DATA_OMICS_PROJECT/DATA"
INPUT_DATA="$PROJECT_DIR/cleaned_data/wakanda_harmonized"
OUT_DIR="$PROJECT_DIR/results/sensitivity"
# crreate output directory if it doesnt exist
mkdir -p "$OUT_DIR"

echo "--- 1. Generating Raw QC Metrics ---"
# We compute missingness, het, and sex check all at once
plink --bfile "$INPUT_DATA" \
      --missing \
      --het \
      --check-sex \
      --hardy \
      --out "$OUT_DIR/raw_qc_metrics"

# ---------------------------------------------------------
# 2. RELATEDNESS CHECK  
# ---------------------------------------------------------
echo "--- Step 2: Preparing Pruned Dataset for Relatedness ---"

# A. Identify Independent SNPs (LD Pruning) on Autosomes
plink --bfile "$INPUT_DATA" \
      --chr 1-22 \
      --indep-pairwise 50 5 0.2 \
      --out "$OUT_DIR/indep" > /dev/null

# B. Create a physical "Mini-Dataset" (Pruned Data)
# This creates a tiny .bed file with only ~30k SNPs.
# It makes the subsequent loop lightning fast and prevents syntax errors.
plink --bfile "$INPUT_DATA" \
      --extract "$OUT_DIR/indep.prune.in" \
      --make-bed \
      --out "$OUT_DIR/pruned_data" > /dev/null

echo "--- Step 3: Running Relatedness Loop on Pruned Data ---"
echo -e "Threshold\tSurvivors" > "$OUT_DIR/sensitivity_relatedness.txt"

REL_THRESHOLDS=("0.05" "0.1" "0.125" "0.15" "0.2" "0.25" "0.3")

for T in "${REL_THRESHOLDS[@]}"
do
    # Run the check on the small 'pruned_data' file
    # We do NOT generate a new bed file (--make-bed), we just check the log.
    plink --bfile "$OUT_DIR/pruned_data" \
          --rel-cutoff "$T" \
          --out "$OUT_DIR/temp_rel_${T}" > /dev/null

    # Extract survivor count
    if [ -f "$OUT_DIR/temp_rel_${T}.log" ]; then
        # Grep for "people remaining"
        REMAINING=$(grep "people remaining" "$OUT_DIR/temp_rel_${T}.log" | awk '{print $1}')
        if [ -z "$REMAINING" ]; then REMAINING="Error"; fi
    else
        REMAINING="Error"
    fi

    echo -e "$T\t$REMAINING" >> "$OUT_DIR/sensitivity_relatedness.txt"
    echo "      Threshold $T: $REMAINING survivors"

    # Clean up temp files
    rm "$OUT_DIR/temp_rel_${T}"* 2>/dev/null
done

# Clean up the mini-dataset
rm "$OUT_DIR/pruned_data"* "$OUT_DIR/indep"*

echo "========================================================"
echo " DONE. Now run the Python plotting script."
echo "========================================================"
