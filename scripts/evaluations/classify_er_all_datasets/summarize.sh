#!/bin/bash
#SBATCH --job-name=gold_summary          # Job name
#SBATCH --output=log/gold_summary_%j.log  # Stdout log (%j = job ID)
#SBATCH --error=log/gold_summary_%j.err   # Stderr log
#SBATCH --time=01:00:00                 # Max run time (hh:mm:ss)
#SBATCH --mem=8G                        # Memory per node
#SBATCH --cpus-per-task=2               # Number of CPUs

# Load modules if necessary
module load python

# Run the Python script
python summarize_gold_data.py
