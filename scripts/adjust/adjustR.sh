#!/bin/bash

# Exit immediately if a command exits with a non-zero status.
set -e

# --- Configuration ---

# Path to the R script that performs the adjustment
ADJUST_SCRIPT="/scripts/adjust/adjust.R"
DATA_DIR="/data"

# Define adjusters to run in PARALLEL (all dataset jobs for an adjuster run at once)
ADJUSTERS_PARALLEL=(
    # "combat"
    # "quantile"
)

# Define adjusters to run SEQUENTIALLY (one dataset job at a time for each adjuster)
ADJUSTERS_SEQUENTIAL=(
    # "seurat_scaling"
    # "seurat_integration"
    # "fastMNN"
    "liger"
    "scanorama"
    "scvi"
)

# Define datasets to be processed by all adjusters
DATASETS=("gse49711" "gse20194" "gse24080")

# Define batch columns for each dataset using an associative array
declare -A BATCH_COLS
BATCH_COLS["gse49711"]="meta_Class"
BATCH_COLS["gse20194"]="meta_batch"
BATCH_COLS["gse24080"]="meta_batch"

# Define target columns to preserve for specific combat and fairadapt adjustments
declare -A TARGET_COLS
TARGET_COLS["gse49711"]="meta_INSS_Stage_Split_1_2 meta_INSS_Stage_Split_2_3 meta_INSS_Stage_Split_3_4"
TARGET_COLS["gse20194"]="meta_er_status meta_her2_status meta_pr_status"
TARGET_COLS["gse24080"]="meta_efs_outcome_label meta_os_outcome_label"


# --- Helper Function ---

# Function to run a single adjustment.
# It constructs and executes the Rscript command.
# Note: The '&' for backgrounding is now handled by the calling loop.
# Usage: run_adjust <adjuster> <dataset> [target_cols_string]
run_adjust() {
    local adjuster=$1
    local dataset=$2
    local target_cols=$3
    
    local input_file="${DATA_DIR}/${dataset}/unadjusted.csv"
    local batch_col="${BATCH_COLS[$dataset]}"
    local output_file
    local c_args=""

    # Determine output filename and -c arguments based on whether target_cols are provided
    if [[ -n "$target_cols" ]]; then
        output_file="${DATA_DIR}/${dataset}/${adjuster}_target.csv"
        c_args="-c ${target_cols}"
    else
        output_file="${DATA_DIR}/${dataset}/${adjuster}.csv"
    fi

    printf " -> Processing %s for %s\n" "$dataset" "$adjuster"
    
    # Execute the R script. The calling loop will decide whether to run in background.
    Rscript "$ADJUST_SCRIPT" "$input_file" "$output_file" -a "$adjuster" -b "$batch_col" $c_args
}


# --- Main Execution ---

# Run adjusters that use specific target columns (Combat, FairAdapt, and Limma) in parallel.
printf "\n\033[0;32mAdjusting data while preserving target columns (Combat, FairAdapt & Limma)\033[0m\n"
# for dataset in "${DATASETS[@]}"; do
#     run_adjust "combat" "$dataset" "${TARGET_COLS[$dataset]}" &
#     run_adjust "fairadapt" "$dataset" "${TARGET_COLS[$dataset]}" &
#     run_adjust "limma" "$dataset" "${TARGET_COLS[$dataset]}" &
# done
# wait

# Run adjusters in the PARALLEL list.
for adjuster in "${ADJUSTERS_PARALLEL[@]}"; do
    printf "\n\033[0;32mAdjusting data with %s (in parallel)\033[0m\n" "$adjuster"
    for dataset in "${DATASETS[@]}"; do
        run_adjust "$adjuster" "$dataset" &
    done
    wait
done

# Run adjusters in the SEQUENTIAL list.
for adjuster in "${ADJUSTERS_SEQUENTIAL[@]}"; do
    printf "\n\033[0;32mAdjusting data with %s (sequentially)\033[0m\n" "$adjuster"
    for dataset in "${DATASETS[@]}"; do
        # Run the job and wait for it to complete before starting the next one.
        run_adjust "$adjuster" "$dataset"
    done
done

printf "\n\033[0;32mAll batch adjustments complete.\033[0m\n"
