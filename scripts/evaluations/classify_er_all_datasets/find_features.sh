#!/bin/bash
#SBATCH --job-name=features
#SBATCH --output=log/features/features_%j.log
#SBATCH --array=0-1049
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=10:00:00

# -----------------------------
# CONFIGURATION
# -----------------------------
# Runs permutation importance on each adjusted CSV file as a job array

# Directory containing CSV results
FEATURE_SCRIPT="classify_feature_importance.py"
DATA_DIR="/grphome/grp_batch_effects/data/adjusted_data"
OUT_DIR="/grphome/grp_batch_effects/data/feature_importance"
CSV_LIST="csv_list.txt"

CSV_FILE=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$CSV_LIST")

if [ -z "$CSV_FILE" ]; then
    echo "No CSV for task ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

echo "Processing $CSV_FILE"

pixi run python "$FEATURE_SCRIPT" \
    --csv "$CSV_FILE" \
    --outdir "$OUT_DIR" \
    --metric roc_auc \
    --n_jobs $SLURM_CPUS_PER_TASK

# Submission line:
# sbatch --array=0-$(($(wc -l < csv_list.txt)-1)) run_features_array.sh