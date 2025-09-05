#!/bin/bash
#
# feature_importance.sh
#
# This script runs feature importance analysis for HistGradientBoostingClassifier
# on the gse_20194_62944/weird_mean.csv dataset
#

# Exit immediately if a command exits with a non-zero status.
set -e

# --- Configuration ---
SCRIPT_PATH="$(dirname "$0")/feature_importance_analysis.py"
DATA_PATH="/data/gold/gse_20194_62944/weird_mean.csv"
OUTPUT_DIR="/outputs/metrics/feature_importance"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

printf "\n\033[0;32mStarting Feature Importance Analysis\033[0m\n"
printf "Data: %s\n" "$DATA_PATH"
printf "Output: %s\n" "$OUTPUT_DIR"

# Check if data file exists
if [ ! -f "$DATA_PATH" ]; then
    printf "\033[0;31mError: Data file not found: %s\033[0m\n" "$DATA_PATH"
    exit 1
fi

# Run the feature importance analysis
printf "\n\033[0;34mRunning feature importance analysis...\033[0m\n"
python3 "$SCRIPT_PATH" \
    --data-path "$DATA_PATH" \
    --output-dir "$OUTPUT_DIR" \
    --create-plots

printf "\n\033[0;32mFeature importance analysis complete!\033[0m\n"
printf "Results saved to: %s\n" "$OUTPUT_DIR"

# List the generated files
printf "\n\033[0;34mGenerated files:\033[0m\n"
ls -la "$OUTPUT_DIR"