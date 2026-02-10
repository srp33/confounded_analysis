#!/bin/bash
#SBATCH --job-name=permute_data
#SBATCH --array=0-4
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=10:00:00
#SBATCH --output=logs/permutation_importance/importance_%A_%a.out

set -euo pipefail

# Set up fixed inputs
IMPORTANCE_SCRIPT="permutation_importance.py"
ADJUST_DIR="/grphome/grp_batch_effects/outputs/metadata_features/labeled_adjusted"
OUTPUT_DIR="/grphome/grp_batch_effects/outputs/metadata_features"
RANDOM_STATE=42
N_REPEATS=3
N_JOBS=${SLURM_CPUS_PER_TASK}

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

ADJ=$(basename "$CSV" | sed -E 's/_([0-9]+)studies_test_.*$//')
MODELS_DIR="/grphome/grp_batch_effects/outputs/metadata_features/classify/${ADJ}/models"

mkdir -p "${OUTPUT_DIR}"

echo "======================================"
echo "SLURM job: $SLURM_JOB_ID"
echo "Task ID: $SLURM_ARRAY_TASK_ID"
echo "CSV file: $CSV"
echo "======================================"

# Run permutation importance
pixi run python "${IMPORTANCE_SCRIPT}" \
    --csv "${CSV}" \
    --models_dir "${MODELS_DIR}" \
    --outdir "${OUTPUT_DIR}" \
    --n_repeats "${N_REPEATS}" \
    --n_jobs "${N_JOBS}" \
    --random_state "${RANDOM_STATE}"

echo "Finished permutation importance ${FILE_IDX}."