#!/bin/bash
#
# reduce.sh
#
# This script automates the process of running dimensionality reduction (PCA, LDA, etc.)
# on a set of batch-corrected data files using a Python script.
#
# For each dataset and for each reduction type, it iterates through all the
# corrected CSV files (e.g., combat.csv, unadjusted.csv) and calls the
# Python script to generate a new CSV file containing the 2D coordinates.
# These output files can then be used by a separate plotting script.

# Exit immediately if a command exits with a non-zero status.
set -e

# --- Configuration ---

# Path to the Python script that performs dimensionality reduction
REDUCTION_SCRIPT="/scripts/metrics/reduce.py"

# Base directory where dataset folders containing corrected CSVs are located
DATA_DIR="/data"

# Base directory where output CSVs with reduced dimensions will be saved
OUTPUT_DIR="/data/reduced_data"

# Define datasets to be processed
# Add or remove dataset names as needed.
DATASETS=(
    "gse49711"
    "gse20194"
    "gse24080"
    # "special_distinct"
)

# Define the types of dimensionality reduction to perform
REDUCTION_TYPES=("pca" "lda" "tsne" "umap")

# --- Dataset-Specific Column Configuration ---

# Define batch columns for each dataset using an associative array
declare -A BATCH_COLS
BATCH_COLS["gse49711"]="meta_Class"
BATCH_COLS["gse20194"]="meta_batch"
BATCH_COLS["gse24080"]="meta_batch"
BATCH_COLS["special_distinct"]="batch"

# Define the "true biological signal" column for each dataset.
declare -A TRUE_COLS
TRUE_COLS["gse49711"]="meta_INSS_Stage_Split_1_2"
TRUE_COLS["gse20194"]="meta_er_status"
TRUE_COLS["gse24080"]="meta_efs_outcome_label"
TRUE_COLS["special_distinct"]="class"


# --- Main Execution ---

printf "\n\033[0;32mStarting Dimensionality Reduction with Python\033[0m\n"

# Loop through each dataset
for dataset in "${DATASETS[@]}"; do
    printf "\n\033[0;34mProcessing dataset: %s\033[0m\n" "$dataset"

    # Define paths and column names for the current dataset
    input_dir="${DATA_DIR}/${dataset}"
    batch_col="${BATCH_COLS[$dataset]}"
    true_col="${TRUE_COLS[$dataset]}"
    output_dir_dataset="${OUTPUT_DIR}/${dataset}"

    # Check if the required columns are defined for the dataset
    if [[ -z "$batch_col" || -z "$true_col" ]]; then
        printf "\033[0;31mWarning: Batch or True column not defined for %s. Skipping.\033[0m\n" "$dataset"
        continue
    fi

    # Create the output directory for the current dataset
    mkdir -p "$output_dir_dataset"

    # Loop through each reduction type for the current dataset
    for reduction_type in "${REDUCTION_TYPES[@]}"; do
        printf " -> Applying %s reduction...\n" "$(echo "$reduction_type" | tr '[:lower:]' '[:upper:]')"

        # Find all corrected CSV files in the input directory
        for input_file in "$input_dir"/*.csv; do
            if [ -f "$input_file" ]; then
                # Extract the method name from the filename (e.g., "combat" from "combat.csv")
                method_name=$(basename "$input_file" .csv)
                printf "    - Processing method: %s\n" "$method_name"

                # Define the output file path
                output_file="${output_dir_dataset}/${method_name}_${reduction_type}.csv"

                # Execute the Python script with file-specific arguments.
                # The script is expected to take one input CSV, run reduction,
                # and save the coordinates to the output CSV.
                python3 "$REDUCTION_SCRIPT" \
                    --input-file "$input_file" \
                    --output-file "$output_file" \
                    --batch-col "$batch_col" \
                    --true-col "$true_col" \
                    --reduction-type "$reduction_type" \
                    --debug
            fi
        done
    done
done

printf "\n\033[0;32mAll dimensionality reduction tasks are complete.\033[0m\n"
