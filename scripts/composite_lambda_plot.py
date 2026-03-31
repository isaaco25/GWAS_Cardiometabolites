import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import re

# Input/Output Files
summary_file = "results/gwas_results/stepwise_analysis/lambda_summary_report.txt"
output_plot = "results/gwas_results/stepwise_analysis/all_plots/Lambda_Trend_Summary.png"

# Read the data
print(f"Reading {summary_file}...")
try:
    df = pd.read_csv(summary_file, sep="\t")
except:
    # Try reading with varying whitespace if tab fails
    df = pd.read_csv(summary_file, delim_whitespace=True)

# Function to extract Step Number (e.g., "step3" -> 3)
def get_step_num(step_name):
    match = re.search(r'step(\d+)', step_name)
    if match:
        return int(match.group(1))
    return 0

# Apply extraction
df['Step_Num'] = df['Step'].apply(get_step_num)

# Sort by Step Number
df = df.sort_values('Step_Num')

# Create the Plot
plt.figure(figsize=(12, 8))
sns.set_style("whitegrid")

# Draw a line for each Trait
# We only plot traits that appear in Step 3 (to avoid clutter if some are missing)
unique_traits = df['Trait'].unique()

print("Plotting trends...")
sns.lineplot(data=df, x='Step_Num', y='Lambda_GC', hue='Trait', marker='o', linewidth=2.5)

# Add a red dashed line at Lambda = 1.0 (The Goal)
plt.axhline(y=1.0, color='r', linestyle='--', linewidth=2, label="Perfect (1.0)")

# Formatting
plt.title("Lambda Inflation Factor by Step (How many PCs do we need?)", fontsize=16)
plt.xlabel("Step Number (Step 1=Base, Step 3=1PC, Step 12=10PCs)", fontsize=12)
plt.ylabel("Lambda (Inflation Factor)", fontsize=12)
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()

# Save
plt.savefig(output_plot, dpi=300)
print(f"Composite Plot saved to: {output_plot}")
