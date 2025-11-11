#!/bin/bash
#SBATCH --job-name=run_classifier
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=logs/classifier_%A_%a.out

# /scripts/evaluations/classify_er_all_datasets/run_classifier.sh

# --- Set up directories ---
CLASSIFIER_SCRIPT="/scripts/evaluations/classify_er_all_datasets/run_classifier.py"
DATA_DIR="/data/adjusted_datasets"
OUTPUT_DIR="/outputs/classify_er_all"

# --- Validate SLURM_ARRAY_TASK_ID ---
if [ -z "$SLURM_ARRAY_TASK_ID" ]; then
    echo "ERROR: SLURM_ARRAY_TASK_ID is not set. Submit this script as a job array." >&2
    exit 1
fi

# --- Folder Logic ---
ADJUSTERS=("gmm" "log_combat" "min_mean" "mnn")
ADJ_FOLDER="${ADJUSTERS[$SLURM_ARRAY_TASK_ID]}"
CSV_DIR="$DATA_DIR/$ADJ_FOLDER"
OUTPUT_CSV="$OUTPUT_DIR/${ADJ_FOLDER}.csv"

mkdir -p "$OUTPUT_DIR"

# Write CSV header
echo "dataset,accuracy,precision,recall,f1" > "$OUTPUT_CSV"

# Check if folder is empty
shopt -s nullglob
CSV_FILES=("$CSV_DIR"/*.csv)

if [ ${#CSV_FILES[@]} -eq 0 ]; then
    echo "No CSV files found in $CSV_DIR"
    exit 1
fi 

# --- Run classifier on each CSV in the folder ---
for CSV_FILE in "${CSV_FILES[@]}"; do
    # Run classifier and capture metrics
    python "$CLASSIFIER_SCRIPT" "$CSV_FILE" "$OUTPUT_DIR"
    STATUS=$?


    if [ $STATUS -eq 0 ]; then
        echo "Classifier finished successfully for $CSV_FILE"

    else
        echo "ERROR: Classifier failed for $CSV_FILE" >&2
    fi
done