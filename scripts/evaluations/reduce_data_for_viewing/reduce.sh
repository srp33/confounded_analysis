#!/bin/bash
#
# reduce.sh
#
# This script automates the process of running dimensionality reduction
# on a set of batch-corrected data files using a Python script.
# It delegates caching to the Python script to avoid re-processing.

# Exit immediately if a command exits with a non-zero status.
set -e

# --- Configuration ---

# Path to the Python script that performs dimensionality reduction
# (This script is expected to handle its own hash-based caching)
REDUCTION_SCRIPT="$(dirname "$0")/reduce.py"

# Base directory where dataset folders are located
DATA_DIR="/data/gold"

# Base directory for output files
OUTPUT_DIR="/data/gold/reduced_data"

# Central directory to store cache files
HASH_DIR="/data/.cache"

# Define datasets to be processed
DATASETS=(
    "gse20194"
    "gse49711"
    "gse24080"
    # "special_distinct"
)

# --- Dataset-Specific Column Configuration ---

declare -A BATCH_COLS
BATCH_COLS["gse49711"]="meta_Sex"
BATCH_COLS["gse20194"]="meta_batch"
BATCH_COLS["gse24080"]="meta_batch"
BATCH_COLS["special_distinct"]="batch"

# TRUE_COLS: Primary biological signal columns used for LDA (Linear Discriminant Analysis)
# LDA requires class labels for supervised dimensionality reduction
declare -A TRUE_COLS
TRUE_COLS["gse49711"]="meta_INSS_Stage_Split_3_4"
TRUE_COLS["gse20194"]="meta_er_status"
TRUE_COLS["gse24080"]="meta_efs_outcome_label"
TRUE_COLS["special_distinct"]="class"

# Additional metadata columns to include in the output
declare -A ADDITIONAL_META_COLS
ADDITIONAL_META_COLS["gse49711"]=""
ADDITIONAL_META_COLS["gse20194"]="meta_her2_status meta_pr_status"
ADDITIONAL_META_COLS["gse24080"]="meta_cytogenetic_abnormality"
ADDITIONAL_META_COLS["special_distinct"]=""


# --- Main Execution ---

printf "\n\033[0;32mStarting Dimensionality Reduction\033[0m\n"

# Create the central cache directory if it doesn't exist
mkdir -p "$HASH_DIR"

# Loop through each dataset
for dataset in "${DATASETS[@]}"; do
    printf "\n\033[0;34mProcessing dataset: %s\033[0m\n" "$dataset"

    # Define paths and column names for the current dataset
    input_dir="${DATA_DIR}/${dataset}"
    batch_col="${BATCH_COLS[$dataset]}"
    true_col="${TRUE_COLS[$dataset]}"
    additional_meta_cols="${ADDITIONAL_META_COLS[$dataset]}"
    output_dir_dataset="${OUTPUT_DIR}/${dataset}"

    # Check if the required columns are defined for the dataset
    if [[ -z "$batch_col" || -z "$true_col" ]]; then
        printf "\033[0;31mWarning: Batch or True column not defined for %s. Skipping.\033[0m\n" "$dataset"
        continue
    fi

    # Create the output directory for the dataset
    mkdir -p "$output_dir_dataset"

    # Find all corrected CSV files in the input directory
    for input_file in "$input_dir"/*.csv; do
        # Process the file only if it exists and is not a transposed file (ending in _t.csv)
        if [ -f "$input_file" ] && [[ "$input_file" != *_t.csv ]]; then
            ( # Start a subshell for each parallel job
                # Build the command with optional additional metadata columns
                cmd_args=(
                    python3 "$REDUCTION_SCRIPT"
                    --input-file "$input_file"
                    --output-dir "$output_dir_dataset"
                    --batch-col "$batch_col"
                    --true-col "$true_col"
                    --hash-dir "$HASH_DIR"
                    # --write-over
                    --debug
                )
                
                # Add additional metadata columns if they exist
                if [[ -n "$additional_meta_cols" ]]; then
                    cmd_args+=(--additional-meta-cols $additional_meta_cols)
                fi
                
                # Execute the command
                "${cmd_args[@]}"
            ) & # Run the subshell in the background
        fi
    done
done

wait # Wait for all background processes to finish

printf "\n\033[0;32mAll dimensionality reduction tasks are complete.\033[0m\n"
