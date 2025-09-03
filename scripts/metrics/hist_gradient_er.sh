#!/bin/bash

# Exit immediately if a command exits with a non-zero status.
set -e

printf "\033[0;32mRunning Classification experiment\033[0m\n"

# --- Configuration ---
# Define paths and parameters used in the script.
script_path="/scripts/metrics/hist_gradient_er_classification.py"
input_base="/data/combined_data/gse20194_gse96058_hiseq"
out_path="/outputs/metrics/hist_gradient_er_results.csv"
summary_path="/outputs/metrics/hist_gradient_er_summary_gse20194_gse96058_hiseq.csv"
cache_dir="/data/.cache/"
# Standardize the number of repeats for all experiments.
n_repeats=2

# --- Filename to Adjustment Mapping ---
# Use an associative array to define the mapping from the input filename
# (without extension) to the adjustment name used in the report.
declare -A adjustment_map
adjustment_map=(
    ["unadjusted"]="Unadjusted"
    ["combat"]="Combat"
    ["quantile"]="Quantile"
    ["min_mean"]="Min-Mean"
    ["fastMNN"]="fastMNN"
    ["npn"]="NPN"
    ["combat_target"]="Combat Target"
    ["simple"]="Simple"
    ["gmm"]="GMM"
    ["gmm_scale_separate"]="GMM Scale Separate"
    ["gmm_npn"]="GMM NPN"
    ["gmm_npn_unit_std"]="GMM NPN Unit Std"
)

# --- Main loop ---
# Find all .csv files in the input directory and iterate over them.
# This avoids the need for a hardcoded execution order.
for input_file in $(find "${input_base}" -name "*.csv"); do
    # Extract the base filename without the path and .csv extension.
    filename_base=$(basename "${input_file}" .csv)

    # Look up the styled adjustment name from the map.
    adjustment="${adjustment_map[$filename_base]}"

    # If a mapping doesn't exist, use the filename base as the adjustment name.
    if [[ -z "${adjustment}" ]]; then
        printf "\n\033[0;33mInfo: No adjustment mapping found for ${filename_base}.csv. Using filename as adjustment name.\033[0m\n"
        adjustment="${filename_base}"
    fi

    printf "\n\033[0;34mRunning for adjustment: ${adjustment} (file: ${filename_base}.csv)\033[0m\n"

    # Execute the python script with the current adjustment's parameters.
    python "${script_path}" \
        --input-data "${input_file}" \
        --summary "${summary_path}" \
        --output "${out_path}" \
        --adjustment "${adjustment}" \
        --cache-dir "${cache_dir}" \
        --n-repeats "${n_repeats}" \
        --force-rerun \

done

printf "\n\033[0;32mAll experiments finished successfully.\033[0m\n"
