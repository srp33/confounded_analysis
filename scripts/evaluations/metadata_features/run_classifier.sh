#!/bin/bash
#SBATCH --job-name=classify_data
#SBATCH --array=0-49 # (5 adjusted files * 10 chunks) - 1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=logs/classify/classifier_%A_%a.out

set -euo pipefail

# Set up fixed inputs
CLASSIFY_SCRIPT="../classify_er_all_datasets/run_classifier.py"
ADJUSTED_DIR="/grphome/grp_batch_effects/outputs/metadata_features/adjusted"
OUTPUT_DIR="/grphome/grp_batch_effects/outputs/metadata_features/classify"
CHUNK_SIZE=10

# Collect all adjusted CSVs
mapfile -t CSV_FILES < <(find "${ADJUSTED_DIR}" -name "*.csv" | sort)

NUM_FILES=${#CSV_FILES[@]}
NUM_CHUNKS=10 # per file

TOTAL_TASKS=$((NUM_FILES * NUM_CHUNKS))

if [[ "$SLURM_ARRAY_TASK_ID" -ge "$TOTAL_TASKS" ]]; then
    echo "Task ${SLURM_ARRAY_TASK_ID} exceeds total tasks ${TOTAL_TASKS}"
    exit 0
fi

# Decode task index
FILE_IDX=$((SLURM_ARRAY_TASK_ID / NUM_CHUNKS))
CHUNK_IDX=$((SLURM_ARRAY_TASK_ID % NUM_CHUNKS))

CSV="${CSV_FILES[$FILE_IDX]}"

mkdir -p "${OUTPUT_DIR}"

echo "======================================"
echo "SLURM job: $SLURM_JOB_ID"
echo "Task ID: $SLURM_ARRAY_TASK_ID"
echo "CSV file: $CSV"
echo "Chunk: $CHUNK_IDX"
echo "======================================"

# Run classifier
pixi run python "${CLASSIFY_SCRIPT}" \
    --csv "${CSV}" \
    --outdir "${OUTPUT_DIR}" \
    --chunk "${CHUNK_IDX}" \
    --chunk-size "${CHUNK_SIZE}"

echo "Finished classifier chunk ${CHUNK_IDX}."