#!/bin/bash

#SBATCH --time=160:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=4G
#SBATCH --array=1-16
#SBATCH -J "simDEpipe"
#SBATCH -o ./DE_large/logs/simDEpipe_%A_%a.out
#SBATCH -e ./DE_large/logs/err_simDEpipe_%A_%a.err

echo "Running simulation DE pipeline"
echo "Job ID: $SLURM_JOB_ID"
echo "Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "Date: $(date)"

# Define parameter arrays
mean_seq=(1 1.5 2 3)
disp_seq=(1 2 3 4)

# Calculate which parameter combination to use based on array task ID
# Array indices: 1-16 for 4x4 parameter grid
task_id=$((SLURM_ARRAY_TASK_ID - 1))  # Convert to 0-based indexing
mean_idx=$((task_id / 4))
disp_idx=$((task_id % 4))

mean_val=${mean_seq[$mean_idx]}
disp_val=${disp_seq[$disp_idx]}

# Create experiment name matching original format
curr_exp_name="largeM$(echo "$mean_val * 10" | bc | cut -d. -f1)D${disp_val}"

echo "Running experiment: $curr_exp_name"
echo "Mean: $mean_val, Dispersion: $disp_val"

# Set up working directory
cd $(dirname "$0")

# Run the R script with parameters
Rscript sim_DEpipe.R $mean_val $disp_val 20

echo "Completed experiment: $curr_exp_name"
echo "End time: $(date)"