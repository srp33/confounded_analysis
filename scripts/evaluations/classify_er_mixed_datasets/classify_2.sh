#!/bin/bash

# /scripts/evaluations/classify_er_mixed_datasets/classify.sh
# Given an adjuster name, create lists of files from dataset directories 
# in /data/combined_data and pass them to HGB_classify_er.py

# Exit immediately if a command exits with a non-zero status
set -e

# --- Configuration ---
COMBINED_DATA_DIR="/data/combined_data"
SINGLE_DATA_DIR="/data/gold"
SCRIPT_PATH="$(dirname "$0")/HGB_classify_er.py"
OUTPUT_DIR="/outputs/metrics"

# Get adjuster from argument
ADJUSTER="${1}"

# --- Main Execution ---
echo "Current working directory: $(pwd)"
printf "\n\033[0;32mClassifying ER status with adjuster: %s\033[0m\n" "$ADJUSTER"

# Find all dataset files for the given adjuster
echo "Searching for datasets in $COMBINED_DATA_DIR with adjuster: $ADJUSTER"
readarray -t COMBINED_FILES < <(find "$COMBINED_DATA_DIR" -mindepth 1  -name "${ADJUSTER}.csv" -type f)

echo "Found ${#COMBINED_FILES[@]} dataset files"
echo "Combined datasets: ${#COMBINED_FILES[@]}"

# echo "Single datasets: ${#SINGLE_FILES[@]}"

# Prepare output file paths
OUTPUT_FILE="${OUTPUT_DIR}/er_classification_${ADJUSTER}.csv"

# Debug: Print file lists
# echo "Combined files:"
# printf '%s\n' "${COMBINED_FILES[@]}"
# echo "Single files:"
# printf '%s\n' "${SINGLE_FILES[@]}"

# Run the Python classification script
echo "Running HGB classification..."
python "$SCRIPT_PATH" \
    --combined-list "${COMBINED_FILES[@]}" \
    --single-list "${SINGLE_FILES[@]}" \
    --output "$OUTPUT_FILE" \
    --adjustment "$ADJUSTER" \
    --classifier "HistGradientBoosting" \
    --prediction-column "meta_er_status" \
    --source-column "meta_source" \
    --n 10