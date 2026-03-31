#!/bin/bash
# diagnose_relatedness_only.sh

# 2. CONFIG
PROJECT_DIR="/home/abegunde/DATA_OMICS_PROJECT/DATA"
INPUT_DATA="$PROJECT_DIR/cleaned_data/wakanda_harmonized"
OUT_DIR="$PROJECT_DIR/results/relatedness_diagnostic"
mkdir -p "$OUT_DIR"

echo "========================================================"
echo " RUNNING ROBUST RELATEDNESS DIAGNOSTIC"
echo "========================================================"

# 3. PREPARE PRUNED MINI-DATASET 
echo "--- Creating Pruned Mini-Dataset ---"
plink --bfile "$INPUT_DATA" \
      --chr 1-22 \
      --maf 0.05 \
      --indep-pairwise 50 5 0.2 \
      --out "$OUT_DIR/indep" > /dev/null

plink --bfile "$INPUT_DATA" \
      --extract "$OUT_DIR/indep.prune.in" \
      --make-bed \
      --out "$OUT_DIR/mini_data" > /dev/null

# 4. RUN THE LOOP
# We test a full spectrum
THRESHOLDS=("0.05" "0.1" "0.125" "0.2" "0.3" "0.4" "0.5" "0.6" "0.7" "0.8" "0.9" "0.95")
OUTPUT_FILE="$OUT_DIR/survivor_counts.txt"

echo -e "Threshold\tSurvivors" > "$OUTPUT_FILE"

echo "--- Testing Thresholds ---"
for T in "${THRESHOLDS[@]}"
do
    # Run cutoff
    plink --bfile "$OUT_DIR/mini_data" \
          --rel-cutoff "$T" \
          --make-bed \
          --out "$OUT_DIR/temp_${T}" > /dev/null
    
    # COUNT SURVIVORS DIRECTLY (Count lines in .fam file)
    # This is 100% accurate and ignores log file text variations.
    if [ -f "$OUT_DIR/temp_${T}.fam" ]; then
        COUNT=$(wc -l < "$OUT_DIR/temp_${T}.fam")
    else
        COUNT="0"
    fi
    
    echo -e "$T\t$COUNT" >> "$OUTPUT_FILE"
    echo "   -> Threshold $T : $COUNT people kept"
    
    # Cleanup temp files
    rm "$OUT_DIR/temp_${T}"* 2>/dev/null
done

# Cleanup mini dataset
rm "$OUT_DIR/mini_data"* "$OUT_DIR/indep"*

echo "========================================================"
echo " DONE. Data saved to: $OUTPUT_FILE"
echo "========================================================"
