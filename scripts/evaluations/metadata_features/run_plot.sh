#!/bin/bash
#SBATCH --job-name=plot_and_select
#SBATCH --array=0-4
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=logs/permutation_importance/importance_%A_%a.out

set -euo pipefail

# Set up fixed inputs
IMPORTANCE_SCRIPT="permutation_importance.py"
ADJUST_DIR="/grphome/grp_batch_effects/outputs/metadata_features/labeled_adjusted"
OUTPUT_DIR="/grphome/grp_batch_effects/outputs/metadata_features"
TOP_K=100
THRESHOLD=0.005

# Collect all adjusted CSVs
mapfile -t CSV_FILES < <(find "${ADJUST_DIR}" -name "*.csv" | sort)

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
pixi run python "${IMPORTANCE_SCRIPT}" \
    --importance_csv "${CSV}" \
    --outdir "${OUTPUT_DIR}" \
    --top_k "${N_REPEATS}" \
    --threshold "${}

echo "Finished permutation importance ${FILE_IDX}."