#!/bin/bash
#SBATCH --job-name=Manhattan_Plots
#SBATCH --mem=16gb             # Give it 16GB RAM (plotting takes memory!)
#SBATCH --output=logs/Plots_%j.log

# 1. Load the R software module (The Cluster needs this first)
module load R

# 2. Run the R script 

Rscript scripts/manhattan_plotsv2.R
