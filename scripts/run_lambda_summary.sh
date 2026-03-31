#!/bin/bash

# ==============================================================================
# CONFIGURATION
# ==============================================================================
# 1. Where are the step folders located?
RESULTS_DIR="/lscratch/GWAS_practice/DATA/results/gwas_results/stepwise_analysis"

# 2. Where is the Python plotting script?
PYTHON_SCRIPT="/lscratch/GWAS_practice/DATA/scripts/lambda_plotv2.py"

# 3. Where should we save the final summary report?
REPORT_FILE="${RESULTS_DIR}/lambda_summary_report.txt"

# ==============================================================================
# EXECUTION
# ==============================================================================
echo "Starting Lambda Analysis..."
echo "Looking for step folders in: $RESULTS_DIR"

# Initialize the Report File with a Header
echo -e "Step\tTrait\tLambda_GC" > "$REPORT_FILE"

# Move to the results directory
cd "$RESULTS_DIR" || { echo "Error: Could not find results directory"; exit 1; }

# Loop through every folder (assuming they are named step1, step2, etc.)
for DIR in step*; do
    if [ -d "$DIR" ]; then
        echo "Processing directory: $DIR"
        
        # Enter the directory
        cd "$DIR"
        
        # Loop through all linear association files
        for FILE in *.assoc.linear; do
            if [ -f "$FILE" ]; then
                # Define output image name (e.g., assoc_bmi.png)
                OUT_IMG="${FILE}.png"
                
                # Run the Python script
                # We capture the output text to grab the Lambda number
                SCRIPT_OUTPUT=$(python "$PYTHON_SCRIPT" "$FILE" "$OUT_IMG")
                
                # Extract the Lambda value using grep and awk
                # It looks for the line "Lambda (GC): 1.0XXX"
                LAMBDA_VAL=$(echo "$SCRIPT_OUTPUT" | grep "Lambda (GC):" | awk '{print $NF}')
                
                # Clean up the trait name for the report
                TRAIT_NAME=$(echo "$FILE" | sed 's/assoc_//' | sed 's/\.assoc\.linear//')
                
                # Print to terminal for progress monitoring
                echo "   -> $TRAIT_NAME: $LAMBDA_VAL"
                
                # Save to the summary report
                echo -e "${DIR}\t${TRAIT_NAME}\t${LAMBDA_VAL}" >> "$REPORT_FILE"
            fi
        done
        
        # Go back up to the main folder
        cd "$RESULTS_DIR"
    fi
done

echo "========================================================"
echo "Analysis Complete!"
echo "Summary Report saved to: $REPORT_FILE"
echo "QQ Plots saved inside each step folder."
echo "========================================================"
