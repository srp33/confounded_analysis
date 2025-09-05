#!/bin/bash
#
# plot_reduced.sh
#
# This script calls the R plotting script for each dataset. The R script
# handles all internal logic, including file discovery and plotting.

set -e

# --- Configuration ---
PLOT_SCRIPT="$(dirname "$0")/plot_reduced.R"
BASE_INPUT_DIR="/data/gold/reduced_data"
OUTPUT_DIR="/outputs/figures/reduced"

# Define datasets to be processed
DATASETS=(
    "gse20194"
    "gse49711"
    "gse24080"
)

# --- Main Execution ---
printf "\n\033[0;32mStarting R Plot Generation\033[0m\n"

# Ensure the output directory exists
mkdir -p "$OUTPUT_DIR"

# Loop through each dataset
for dataset in "${DATASETS[@]}"; do
    printf "\n\033[0;34mProcessing dataset: %s\033[0m\n" "$dataset"

    # Define paths for the current dataset
    input_dir="${BASE_INPUT_DIR}/${dataset}"
    output_pdf="${OUTPUT_DIR}/${dataset}.pdf"

    # Execute the R script ONCE per dataset.
    # The R script will find all relevant files and generate one PDF.
    Rscript "$PLOT_SCRIPT" \
        -i "$input_dir" \
        -o "$output_pdf"
done

printf "\n\033[0;32mAll R plot generation tasks are complete.\033[0m\n"