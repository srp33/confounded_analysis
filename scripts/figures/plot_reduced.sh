#!/bin/bash
#
# run_plots_r.sh
#
# This script automates the generation of dimensionality reduction plots
# using the R plotting script (plot_reduced.R).
#
# It iterates through a predefined list of datasets and plot types (PCA, LDA, etc.),
# calling the R script for each combination. The R script itself is
# internally parallelized to process all correction methods for a given dataset.

# Exit immediately if a command exits with a non-zero status.
set -e

# --- Configuration ---

# Path to the R script that generates the plots
PLOT_SCRIPT="/scripts/figures/plot_reduced.R"

# Base directory where dataset folders are located
DATA_DIR="/data"

# Base directory where output plots will be saved
OUTPUT_DIR="/outputs/figures/reduced/"

# Define datasets to be processed
# Add or remove dataset names as needed.
DATASETS=(
    "gse49711"
    "gse20194"
    "gse24080"
    # "special_distinct"
)

# Define the types of dimensionality reduction plots to generate
PLOT_TYPES=("pca" "lda" "tsne" "umap")

# --- Dataset-Specific Column Configuration ---

# Define batch columns for each dataset using an associative array
declare -A BATCH_COLS
BATCH_COLS["gse49711"]="meta_Class"
BATCH_COLS["gse20194"]="meta_batch"
BATCH_COLS["gse24080"]="meta_batch"
BATCH_COLS["special_distinct"]="batch"

# Define the "true biological signal" column for each dataset.
# This is used for the 'shape' aesthetic in the plots.
declare -A TRUE_COLS
TRUE_COLS["gse49711"]="meta_INSS_Stage_Split_1_2"
TRUE_COLS["gse20194"]="meta_er_status"
TRUE_COLS["gse24080"]="meta_efs_outcome_label"
TRUE_COLS["special_distinct"]="class"


# --- Main Execution ---

printf "\n\033[0;32mStarting R Plot Generation\033[0m\n"

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

    # Loop through each plot type for the current dataset
    for plot_type in "${PLOT_TYPES[@]}"; do
        printf " -> Generating %s plot...\n" "$(echo "$plot_type" | tr '[:lower:]' '[:upper:]')"

        # Create the output directory if it doesn't exist
        mkdir -p "$output_dir_dataset"

        # Execute the R script with all required arguments.
        # The script is expected to find all '*.csv' files in the input_dir,
        # run the reduction, and save a single PDF with all the plots.
        Rscript "$PLOT_SCRIPT" \
            -i "$input_dir" \
            -o "$output_dir_dataset" \
            -b "$batch_col" \
            -t "$true_col" \
            -d "$dataset" \
            -p "$plot_type" \
            --debug
    done
done

printf "\n\033[0;32mAll R plot generation tasks are complete.\033[0m\n"
