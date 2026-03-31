import pandas as pd
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
import os
import glob
import re

# ==============================================================================
# CONFIGURATION
# ==============================================================================
RESULTS_DIR = "/lscratch/GWAS_practice/DATA/results/gwas_results/stepwise_analysis"
OUT_PLOT = "/lscratch/GWAS_practice/DATA/results/gwas_results/lambda_decay_plot.png"
OUT_TABLE = "/lscratch/GWAS_practice/DATA/results/gwas_results/lambda_stats.csv"

# ==============================================================================
# FUNCTIONS
# ==============================================================================
def calculate_lambda(file_path):
    """Calculates Lambda GC from a PLINK .assoc.linear file"""
    try:
        # Read file (space separated), handle errors
        df = pd.read_csv(file_path, sep=r"\s+", engine='python')
        
        # Check if 'P' column exists
        if "P" not in df.columns:
            return None

        # Drop NAs and non-numeric P-values
        df = df.dropna(subset=["P"])
        df = df[pd.to_numeric(df['P'], errors='coerce').notnull()]
        
        # If no valid P-values remain, return None (Skip this file)
        if len(df) == 0:
            return None
        
        # Convert P-values to Chi-Squared statistics
        chi2_stats = stats.chi2.ppf(1 - df["P"], 1)
        
        # Calculate Lambda
        # Median observed / Median expected (0.4549 is median of chi2 with 1 df)
        lambda_gc = np.median(chi2_stats) / 0.454936421
        return lambda_gc
    except Exception as e:
        # Silently fail for bad files
        return None

# ==============================================================================
# MAIN LOGIC
# ==============================================================================
print("--- Starting Lambda Analysis (Robust Mode) ---")

data = []

# 1. Find all Step folders
step_folders = glob.glob(os.path.join(RESULTS_DIR, "step*"))
# Sort by step number
step_folders.sort(key=lambda x: int(re.search(r'step(\d+)', x).group(1)))

print(f"Found {len(step_folders)} step folders.")

# 2. Iterate through steps
for folder in step_folders:
    folder_name = os.path.basename(folder)
    
    match = re.search(r'step(\d+)_added_(.+)', folder_name)
    if not match:
        continue
        
    step_num = int(match.group(1))
    added_covar = match.group(2)
    
    print(f" -> Analyzing Step {step_num}: +{added_covar}")

    # Find all association files in this folder
    assoc_files = glob.glob(os.path.join(folder, "*.assoc.linear"))
    
    for f in assoc_files:
        basename = os.path.basename(f)
        trait_name = basename.replace("assoc_", "").replace(".assoc.linear", "")
        
        lam = calculate_lambda(f)
        
        if lam is not None:
            data.append({
                "Step": step_num,
                "Covariate": added_covar,
                "Trait": trait_name,
                "Lambda": lam
            })

# 3. Create DataFrame
if len(data) == 0:
    print("CRITICAL ERROR: No valid Lambda values could be calculated.")
    exit()

df_results = pd.DataFrame(data)

# Save raw stats to CSV
df_results.sort_values(by=["Trait", "Step"]).to_csv(OUT_TABLE, index=False)
print(f"\nStats saved to: {OUT_TABLE}")

# 4. Plotting
print("Generating Plot...")
plt.figure(figsize=(12, 8))

# Get list of unique traits
traits = df_results["Trait"].unique()

# Plot a line for each trait
for trait in traits:
    subset = df_results[df_results["Trait"] == trait].sort_values("Step")
    plt.plot(subset["Step"], subset["Lambda"], marker='o', label=trait, alpha=0.7)

# Formatting
plt.axhline(y=1.0, color='r', linestyle='--', linewidth=2, label="Ideal (1.0)")
plt.xlabel("Step Number (Cumulative Covariates Added)", fontsize=12)
plt.ylabel("Genomic Inflation Factor (Lambda GC)", fontsize=12)
plt.title("GWAS Inflation Decay: Effect of Adding Covariates", fontsize=14)
plt.grid(True, linestyle='--', alpha=0.5)

# Set X-ticks
step_map = df_results[["Step", "Covariate"]].drop_duplicates().sort_values("Step")
plt.xticks(step_map["Step"], step_map["Covariate"], rotation=45)

plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize='small')
plt.tight_layout()

# Save Plot
plt.savefig(OUT_PLOT, dpi=300)
print(f"Plot saved to: {OUT_PLOT}")
print("Done.")
