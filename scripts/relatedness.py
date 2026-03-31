import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

# CONFIG
file_path = "/home/abegunde/DATA_OMICS_PROJECT/DATA/results/relatedness_diagnostic/survivor_counts.txt"
output_img = "plot_relatedness_final.png"

sns.set_style("whitegrid")

try:
    if not os.path.exists(file_path):
        print(f"Error: File {file_path} not found. Run the bash script first!")
        exit(1)

    # Load Data
    df = pd.read_csv(file_path, sep="\t")
    
    # Plot
    plt.figure(figsize=(12, 7))
    
    # The Main Curve
    plt.plot(df['Threshold'], df['Survivors'], marker='o', linestyle='-', color='navy', linewidth=2.5, label='Survivors')
    
    # Add Data Labels
    for i, txt in enumerate(df['Survivors']):
        plt.annotate(str(txt), (df['Threshold'][i], df['Survivors'][i]), 
                     xytext=(0, 10), textcoords='offset points', ha='center', fontsize=10, fontweight='bold')

    # Add Reference Lines
    plt.axvline(x=0.125, color='red', linestyle='--', linewidth=2, label='0.125 (Remove Cousins)')
    plt.axvline(x=0.9, color='green', linestyle='--', linewidth=2, label='0.9 (Remove Duplicates Only)')
    
    # Add Shaded Regions for Context
    plt.axvspan(0, 0.125, color='red', alpha=0.1)
    plt.axvspan(0.8, 1.0, color='green', alpha=0.1)

    # Labels and Title
    plt.title("Impact of Relatedness Threshold on Sample Size", fontsize=16)
    plt.xlabel("Relatedness Threshold (PI_HAT)\n(Lower = Stricter)", fontsize=12)
    plt.ylabel("Number of Individuals Kept", fontsize=12)
    plt.legend(loc='lower right')
    
    plt.tight_layout()
    plt.savefig(output_img)
    print(f"Plot saved successfully as: {output_img}")

except Exception as e:
    print(f"An error occurred: {e}")
