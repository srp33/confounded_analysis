#!/bin/bash
#
# classify.sh
#

# Exit immediately if a command exits with a non-zero status.
set -e

# Set PYTHONPATH for Python imports
export PYTHONPATH="/scripts:$PYTHONPATH"

# --- Configuration ---
SCRIPT_PATH="$(dirname "$0")/classify.py"
BASE_DATA_DIR="/data/gold"
OUTPUT_DIR="/outputs/metrics"
# Directory to store the cache files with MD5 hashes
HASH_DIR="/data/.cache"

# Create output and cache directories if they don't exist
mkdir -p "$OUTPUT_DIR"
mkdir -p "$HASH_DIR"

# --- Task Definitions ---

# Define datasets and their corresponding true columns for global prediction
declare -A TASKS
TASKS["gse20194"]="meta_er_status meta_her2_status meta_pr_status"
TASKS["gse24080"]="meta_cytogenetic_abnormality"
TASKS["gse49711"]="meta_INSS_Stage_Split_3_4"

TASKS["gse_20194_62944"]="meta_er_status"

TASKS["2_dims_no_bio_no_batch"]="meta_bio"
TASKS["2_dims_no_bio_yes_batch"]="meta_bio"
TASKS["2_dims_yes_bio_no_batch"]="meta_bio"
TASKS["2_dims_yes_bio_yes_batch"]="meta_bio"

TASKS["400_dims_no_bio_no_batch"]="meta_bio"
TASKS["400_dims_no_bio_yes_batch"]="meta_bio"
TASKS["400_dims_yes_bio_no_batch"]="meta_bio"
TASKS["400_dims_yes_bio_yes_batch"]="meta_bio"

TASKS["1000_dims_no_bio_no_batch"]="meta_bio"
TASKS["1000_dims_no_bio_yes_batch"]="meta_bio"
TASKS["1000_dims_yes_bio_no_batch"]="meta_bio"
TASKS["1000_dims_yes_bio_yes_batch"]="meta_bio"

TASKS["structured_synthetic"]="meta_bio"

# Define batch columns for each dataset
declare -A BATCH_COLS
BATCH_COLS["gse20194"]="meta_batch"
BATCH_COLS["gse24080"]="meta_batch"
BATCH_COLS["gse49711"]="meta_Sex"

BATCH_COLS["gse_20194_62944"]="meta_source"

BATCH_COLS["2_dims_no_bio_no_batch"]="meta_batch"
BATCH_COLS["2_dims_no_bio_yes_batch"]="meta_batch"
BATCH_COLS["2_dims_yes_bio_no_batch"]="meta_batch"
BATCH_COLS["2_dims_yes_bio_yes_batch"]="meta_batch"

BATCH_COLS["400_dims_no_bio_no_batch"]="meta_batch"
BATCH_COLS["400_dims_no_bio_yes_batch"]="meta_batch"
BATCH_COLS["400_dims_yes_bio_no_batch"]="meta_batch"
BATCH_COLS["400_dims_yes_bio_yes_batch"]="meta_batch"

BATCH_COLS["1000_dims_no_bio_no_batch"]="meta_batch"
BATCH_COLS["1000_dims_no_bio_yes_batch"]="meta_batch"
BATCH_COLS["1000_dims_yes_bio_no_batch"]="meta_batch"
BATCH_COLS["1000_dims_yes_bio_yes_batch"]="meta_batch"

BATCH_COLS["structured_synthetic"]="meta_batch"


# --- Main Execution ---

printf "\n\033[0;32mStarting Classification Metrics Calculation\033[0m\n"

# == Part 1: Global Predictions (No conditional column) ==
printf "\n\033[0;34mRunning Global Predictions...\033[0m\n"
for dataset in "${!TASKS[@]}"; do
    input_dir="${BASE_DATA_DIR}/${dataset}"
    
    # --- Predict biological true columns ---
    output_path_true="${OUTPUT_DIR}/global_true_classification.csv"
    for p_col in ${TASKS[$dataset]}; do
        python3 "${SCRIPT_PATH}" \
            -i "${input_dir}" \
            -o "${output_path_true}" \
            -p "${p_col}" \
            --hash-dir "${HASH_DIR}"
            # --write-over
    done

    # --- Predict batch column ---
    output_path_batch="${OUTPUT_DIR}/global_batch_classification.csv"
    batch_col="${BATCH_COLS[$dataset]}"
    python3 "${SCRIPT_PATH}" \
        -i "${input_dir}" \
        -o "${output_path_batch}" \
        -p "${batch_col}" \
        --hash-dir "${HASH_DIR}" 
        # --write-over
done

# == Part 2: Conditional Predictions ==
printf "\n\033[0;34mRunning Conditional Predictions...\033[0m\n"
output_path_true="${OUTPUT_DIR}/true_classification.csv"
output_path_batch="${OUTPUT_DIR}/batch_classification.csv"

for dataset in "${!TASKS[@]}"; do
    input_dir="${BASE_DATA_DIR}/${dataset}"
    batch_col="${BATCH_COLS[$dataset]}"

    for p_col in ${TASKS[$dataset]}; do
        # --- Predict true column (P) conditional on batch column (C) ---
        python3 "${SCRIPT_PATH}" -i "${input_dir}" -o "${output_path_true}" -p "${p_col}" -c "${batch_col}" --hash-dir "${HASH_DIR}"

        # --- Predict batch column (P) conditional on true column (C) ---
        python3 "${SCRIPT_PATH}" -i "${input_dir}" -o "${output_path_batch}" -p "${batch_col}" -c "${p_col}" --hash-dir "${HASH_DIR}"
    done
done


printf "\n\033[0;32mAll classification tasks are complete.\033[0m\n"