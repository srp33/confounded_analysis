#!/bin/bash
#SBATCH --job-name=run_classifier
#SBATCH --array=0-3 # One job per adjustser folder
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=logs/classifier_%A_%a.out

# --- Set up directories ---
ANALYSIS_DIR=~/confounded_analysis
CLASSIFIER_SCRIPT=scripts/evaluations/classify_er_all_datasets/run_classifier.py
ADJUSTED_DIR=$ANALYSIS_DIR/grp_batch_effects/data/adjusted_datasets
METRICS_DIR="$ANALYSIS_DIR/grp_batch_effects/data/metrics"

# --- Validate SLURM_ARRAY_TASK_ID ---
if [ -z "$SLURM_ARRAY_TASK_ID" ]; then
    echo "ERROR: SLURM_ARRAY_TASK_ID is not set. Submit this script as a job array." >&2
    exit 1
fi

# --- Folder Logic ---
ADJUSTERS=("gmm" "log_combat" "min_mean" "mnn")
ADJ_FOLDER="${ADJUSTERS[$SLURM_ARRAY_TASK_ID]}"
CSV_DIR="$ADJUSTED_DIR/$ADJ_FOLDER"
OUTPUT_CSV="$METRICS_DIR/${ADJ_FOLDER}.csv"

# Create metrics directory if it doesn't exist
mkdir -p "$METRICS_DIR"

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
    # Original host CSV
    HOST_CSV="$CSV_FILE"

    # Translate host path to container path
    # Assuming /home/aw998/confounded_analysis/grp_batch_effects/data -> /data
    CONTAINER_CSV="${HOST_CSV/\/home\/aw998\/confounded_analysis\/grp_batch_effects\/data/\/data}"
    cd $ANALYSIS_DIR

    echo "[$(date)] Running classifier on host path: $HOST_CSV, container path: $CONTAINER_CSV"

    # Run classifier and capture metrics
    METRICS=$(bash "$ANALYSIS_DIR/run_in_apptainer.sh" "$CLASSIFIER_SCRIPT" "$CONTAINER_CSV")
    STATUS=$?

    if [ $STATUS -eq 0 ]; then
        echo "Classifier finished successfully for $CONTAINER_CSV"
        echo "$(basename $CSV_FILE),$METRICS" >> "$OUTPUT_CSV"

    else
        echo "ERROR: Classifier failed for $CSV_FILE" >&2
    fi
done