import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
import sys
import os

# CONFIGURATION
if len(sys.argv) < 3:
    print("Usage: python plot_qq_lambda.py <input_file> <output_image>")
    sys.exit(1)

input_file = sys.argv[1]
output_file = sys.argv[2]
print(f"Processing: {os.path.basename(input_file)}")

# 1. READ DATA
try:
    df = pd.read_csv(input_file, delim_whitespace=True)
except:
    df = pd.read_csv(input_file, sep='\t')

# 2. FILTER FOR ADDITIVE TEST (Critical for .assoc.linear)
if 'TEST' in df.columns:
    initial_count = len(df)
    # Keep only the genetic effect (ADD)
    df = df[df['TEST'] == 'ADD']
    print(f"   -> Filtered for 'ADD' test: {len(df)} SNPs remain (dropped {initial_count - len(df)} covariate rows)")

# 3. FIND P-VALUE COLUMN
possible_cols = ['P', 'p', 'P_VALUE', 'P.value']
p_col = next((c for c in possible_cols if c in df.columns), None)

if p_col is None:
    print(f"ERROR: No P-value column found in {input_file}")
    sys.exit(1)

# 4. CALCULATE LAMBDA
p_values = df[p_col].dropna()
p_values = p_values[ (p_values > 0) & (p_values <= 1) ]

if len(p_values) == 0:
    print("   -> Error: No valid P-values found.")
    sys.exit(1)

chi2_obs = stats.chi2.ppf(1 - p_values, df=1)
lambda_gc = np.median(chi2_obs) / stats.chi2.ppf(0.5, df=1)
print(f"   -> Lambda (GC): {lambda_gc:.5f}")

# 5. GENERATE PLOT
p_sorted = np.sort(p_values)
n = len(p_sorted)
expected = -np.log10(np.arange(1, n + 1) / (n + 1))
observed = -np.log10(p_sorted)

plt.figure(figsize=(6, 6))
plt.scatter(expected, observed, c='#2c3e50', s=5, alpha=0.6)
plt.plot([0, max(expected)], [0, max(expected)], 'r--', lw=2) # Diagonal

# Add info to plot
pheno_name = os.path.basename(input_file).replace("assoc_", "").replace(".assoc.linear", "")
plt.xlabel("Expected -log10(P)")
plt.ylabel("Observed -log10(P)")
plt.title(f"QQ Plot: {pheno_name}\n$\lambda_{{GC}} = {lambda_gc:.4f}$")
plt.grid(True, alpha=0.3)
plt.tight_layout()

plt.savefig(output_file, dpi=150)
print(f"   -> Plot saved to: {output_file}\n")
