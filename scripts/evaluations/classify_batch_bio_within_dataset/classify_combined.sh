#!/bin/bash
#
# classify_combined.sh - Classification for combined datasets
#

# Exit immediately if a command exits with a non-zero status.
set -e

# --- Configuration ---
SCRIPT_PATH="/scripts/metrics/classify.py"
BASE_DATA_DIR="/data/paired_datasets"
OUTPUT_DIR="/outputs/metrics"
# Directory to store the cache files with MD5 hashes
HASH_DIR="/data/.cache/classify_hashes"

# Create output and cache directories if they don't exist
mkdir -p "$OUTPUT_DIR"
mkdir -p "$HASH_DIR"

# --- Main Execution ---

printf "\n\033[0;32mStarting Classification Metrics Calculation for Combined Data\033[0m\n"

# Get all combined datasets
DATASETS=($(ls -1 "$BASE_DATA_DIR"))

# Classification target and batch column
PREDICT_COL="meta_er_status"
BATCH_COL="meta_source"

# == Part 1: Global Predictions (No conditional column) ==
printf "\n\033[0;34mRunning Global Predictions for Combined Data...\033[0m\n"

for dataset in "${DATASETS[@]}"; do
    input_dir="${BASE_DATA_DIR}/${dataset}"
    
    # Skip if not a directory
    if [ ! -d "$input_dir" ]; then
        continue
    fi
    
    printf "Processing dataset: %s\n" "$dataset"
    
    # --- Predict biological true column ---
    output_path_true="${OUTPUT_DIR}/combined_global_true_classification.csv"
    python3 "${SCRIPT_PATH}" \
        -i "${input_dir}" \
        -o "${output_path_true}" \
        -p "${PREDICT_COL}" \
        --hash-dir "${HASH_DIR}"

    # --- Predict batch column ---
    output_path_batch="${OUTPUT_DIR}/combined_global_batch_classification.csv"
    python3 "${SCRIPT_PATH}" \
        -i "${input_dir}" \
        -o "${output_path_batch}" \
        -p "${BATCH_COL}" \
        --hash-dir "${HASH_DIR}"
done

# == Part 2: Conditional Predictions ==
printf "\n\033[0;34mRunning Conditional Predictions for Combined Data...\033[0m\n"
output_path_true="${OUTPUT_DIR}/combined_true_classification.csv"
output_path_batch="${OUTPUT_DIR}/combined_batch_classification.csv"

for dataset in "${DATASETS[@]}"; do
    input_dir="${BASE_DATA_DIR}/${dataset}"
    
    # Skip if not a directory
    if [ ! -d "$input_dir" ]; then
        continue
    fi
    
    printf "Processing conditional predictions for dataset: %s\n" "$dataset"

    # --- Predict true column (P) conditional on batch column (C) ---
    python3 "${SCRIPT_PATH}" -i "${input_dir}" -o "${output_path_true}" -p "${PREDICT_COL}" -c "${BATCH_COL}" --hash-dir "${HASH_DIR}"

    # --- Predict batch column (P) conditional on true column (C) ---
    python3 "${SCRIPT_PATH}" -i "${input_dir}" -o "${output_path_batch}" -p "${BATCH_COL}" -c "${PREDICT_COL}" --hash-dir "${HASH_DIR}"
done

printf "\n\033[0;32mAll combined data classification tasks are complete.\033[0m\n"