#!/bin/bash
#SBATCH --job-name=adjust_data
#SBATCH --array=1-65%10 # 5 adjusters * 13 subsets = 65 jobs
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --output=logs/adjuster_%A_%a.out

# Set up directory and adjusters
ANALYSIS_DIR=~/confounded_analysis
ADJUSTERS=("mnn" "min_mean" "log_combat" "gmm" "log_transformed")

# Compute adjuster and subset from SLURM_ARRAY_TASK_ID
ADJ_IDX=$(( ($SLURM_ARRAY_TASK_ID - 1) / 13 ))  # 13 subsets per adjuster
SUBSET_IDX=$(( ($SLURM_ARRAY_TASK_ID - 1) % 13 + 2 ))  # subsets 2-14

ADJ=${ADJUSTERS[$ADJ_IDX]}

echo "Requested adjuster: $ADJ on subset: $SUBSET_IDX"

# # Skip all adjusters except log_transformed
# if [ "$ADJ" != "log_transformed" ]; then
#     echo "Skipping adjuster $ADJ (already completed previously)."
#     exit 0
# fi

echo "Running adjuster: $ADJ on subset: $SUBSET_IDX"

cd $ANALYSIS_DIR
bash $ANALYSIS_DIR/run_in_apptainer.sh /scripts/evaluations/classify_er_all_datasets/run_scaling_experiment.R $ADJ $SUBSET_IDX


if [ $? -eq 0 ]; then
    echo "Adjuster $ADJ on subset $SUBSET_IDX finished successfully!"
else
    echo "ERROR: Adjuster $ADJ on subset $SUBSET_IDX failed!" >&2
    exit 1
fi