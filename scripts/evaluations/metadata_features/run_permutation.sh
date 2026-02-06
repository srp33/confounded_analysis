#!/bin/bash
#SBATCH --job-name=classify_data
#SBATCH --array=0-4
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=logs/classify/classifier_%A_%a.out

set -euo pipefail

# Set up fixed inputs
CLASSIFY_SCRIPT="train_and_metrics.py"
ADJUSTED_DIR="/grphome/grp_batch_effects/outputs/metadata_features/labeled_adjusted"
OUTPUT_DIR="/grphome/grp_batch_effects/outputs/metadata_features/classify"
RANDOM_STATE=42
# Collect all adjusted CSVs
mapfile -t CSV_FILES < <(find "${ADJUSTED_DIR}" -name "*.csv" | sort)

NUM_FILES=${#CSV_FILES[@]}

TOTAL_TASKS=$((NUM_FILES))

if [[ "$SLURM_ARRAY_TASK_ID" -ge "$TOTAL_TASKS" ]]; then
    echo "Task ${SLURM_ARRAY_TASK_ID} exceeds total tasks ${TOTAL_TASKS}"
    exit 0
fi

# Decode task index
FILE_IDX=$((SLURM_ARRAY_TASK_ID))

CSV="${CSV_FILES[$FILE_IDX]}"

mkdir -p "${OUTPUT_DIR}"

echo "======================================"
echo "SLURM job: $SLURM_JOB_ID"
echo "Task ID: $SLURM_ARRAY_TASK_ID"
echo "CSV file: $CSV"
echo "======================================"

# Run classifier
pixi run python "${CLASSIFY_SCRIPT}" \
    --csv "${CSV}" \
    --outdir "${OUTPUT_DIR}" \
    --random_state "${RANDOM_STATE}"

echo "Finished classifier ${FILE_IDX}."