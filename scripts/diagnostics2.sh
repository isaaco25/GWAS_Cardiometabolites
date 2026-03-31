#!/bin/bash
# iterative_diagnostics2.sh

# CONFIGURATION
PROJECT_DIR="/home/abegunde/DATA_OMICS_PROJECT/DATA"
INPUT_DATA="$PROJECT_DIR/cleaned_data/wakanda_harmonized"
OUT_DIR="$PROJECT_DIR/results/sensitivity"

mkdir -p "$OUT_DIR"

echo "========================================================"
echo " GENERATING METRICS FOR RELATEDNESS, MAF, AND HWE"
echo "========================================================"

# ---------------------------------------------------------
# 1. MAF & HWE CALCULATIONS
# ---------------------------------------------------------
echo "--- Step 1: Calculating Allele Frequencies & HWE ---"
# We perform this on the full harmonized dataset
plink --bfile "$INPUT_DATA" \
      --freq \
      --hardy \
      --out "$OUT_DIR/raw_qc_metrics"

echo "   -> MAF file created: $OUT_DIR/raw_qc_metrics.frq"
echo "   -> HWE file created: $OUT_DIR/raw_qc_metrics.hwe"

# ---------------------------------------------------------
# 2. RELATEDNESS CHECK (Iterative)
# ---------------------------------------------------------
echo "--- Step 2: Running Iterative Relatedness Check ---"

# A. Create a physical "Mini-Dataset" (Pruned Data) for speed
# We use Autosomes only (1-22) to avoid X-chromosome bias
plink --bfile "$INPUT_DATA" \
      --chr 1-22 \
      --indep-pairwise 50 5 0.2 \
      --out "$OUT_DIR/indep" > /dev/null

plink --bfile "$INPUT_DATA" \
      --extract "$OUT_DIR/indep.prune.in" \
      --make-bed \
      --out "$OUT_DIR/pruned_data" > /dev/null

# B. Loop through thresholds
echo -e "Threshold\tSurvivors" > "$OUT_DIR/sensitivity_relatedness.txt"
REL_THRESHOLDS=("0.05" "0.1" "0.125" "0.2" "0.25" "0.3")

for T in "${REL_THRESHOLDS[@]}"
do
    # Run the check on the small 'pruned_data' file
    plink --bfile "$OUT_DIR/pruned_data" \
          --rel-cutoff "$T" \
          --out "$OUT_DIR/temp_rel_${T}" > /dev/null

    # Extract survivor count
    if [ -f "$OUT_DIR/temp_rel_${T}.log" ]; then
        REMAINING=$(grep "people remaining" "$OUT_DIR/temp_rel_${T}.log" | awk '{print $1}')
        if [ -z "$REMAINING" ]; then REMAINING="Error"; fi
    else
        REMAINING="Error"
    fi

    echo -e "$T\t$REMAINING" >> "$OUT_DIR/sensitivity_relatedness.txt"
    echo "   -> Tested Threshold $T: $REMAINING survivors"

    # Clean up temp files
    rm "$OUT_DIR/temp_rel_${T}"* 2>/dev/null
done

# Clean up intermediate files
rm "$OUT_DIR/pruned_data"* "$OUT_DIR/indep"*

echo "========================================================"
echo " DONE."
echo "========================================================"
