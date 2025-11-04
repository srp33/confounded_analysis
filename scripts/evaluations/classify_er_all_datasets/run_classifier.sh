#!/bin/bash
#SBATCH --job-name=run_classifier
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=logs/classifier_%A_%a.out
#SBATCH --array=1-1%10  # Placeholder, will be overridden by number of files

module load apptainer

CONTAINER=~/confounded_analysis/apptainer/remove-batch-effects.sif
CLASSIFIER_SCRIPT=scripts/evaluations/classify_er_all_datasets/run_classifier.py
ADJUSTED_DIR=grp_batch_effects/data/adjusted_datasets

# --- Detect all CSVs ---
FILES=($(find "$ADJUSTED_DIR" -name "*.csv" | sort))
NUM_FILES=${#FILES[@]}

# --- Update SLURM array dynamically ---
if [ $SLURM_ARRAY_TASK_ID -gt $NUM_FILES ]; then
    echo "No file for this array task ID ($SLURM_ARRAY_TASK_ID), exiting."
    exit 0
fi

CSV_FILE="${FILES[$SLURM_ARRAY_TASK_ID-1]}"

echo "Running classifier on: $CSV_FILE"

apptainer exec "$CONTAINER" python3 "$CLASSIFIER_SCRIPT" --input "$CSV_FILE"

if [ $? -eq 0 ]; then
    echo "Classifier finished successfully for $CSV_FILE"
else
    echo "ERROR: Classifier failed for $CSV_FILE" >&2
    exit 1
fi
