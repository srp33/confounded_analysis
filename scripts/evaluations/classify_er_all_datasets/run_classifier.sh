#!/bin/bash
#SBATCH --job-name=run_classifier
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=10:00:00
#SBATCH --output=logs/classifier_%A_%a.out
#SBATCH --requeue  # allow job to be requeued if killed

# --- Directories ---
CLASSIFIER_SCRIPT="/scripts/evaluations/classify_er_all_datasets/run_classifier.py"
DATA_DIR="/data/adjusted_datasets"
OUTPUT_DIR="/outputs/classify_er_all"

# --- Validate SLURM_ARRAY_TASK_ID ---
if [ -z "$SLURM_ARRAY_TASK_ID" ]; then
    echo "ERROR: SLURM_ARRAY_TASK_ID is not set. Submit this script as a job array." >&2
    exit 1
fi

# --- Adjuster logic ---
ADJUSTERS=("gmm" "log_combat" "min_mean" "mnn")
ADJ_FOLDER="${ADJUSTERS[$SLURM_ARRAY_TASK_ID]}"
CSV_DIR="$DATA_DIR/$ADJ_FOLDER"
OUTPUT_CSV="$OUTPUT_DIR/${ADJ_FOLDER}.csv"

mkdir -p "$OUTPUT_DIR"

# --- CSV header if file does not exist ---
if [ ! -f "$OUTPUT_CSV" ]; then
    echo "dataset,accuracy,precision,recall,f1" > "$OUTPUT_CSV"
fi

# --- Get list of CSVs ---
shopt -s nullglob
CSV_FILES=("$CSV_DIR"/*.csv)

if [ ${#CSV_FILES[@]} -eq 0 ]; then
    echo "No CSV files found in $CSV_DIR"
    exit 1
fi 

# --- Skip datasets that are already processed ---
PROCESSED=$(tail -n +2 "$OUTPUT_CSV" | cut -d',' -f1)  # dataset column
TO_PROCESS=()
for f in "${CSV_FILES[@]}"; do
    fname=$(basename "$f")
    if ! grep -qx "$fname" <<< "$PROCESSED"; then
        TO_PROCESS+=("$f")
    fi
done

if [ ${#TO_PROCESS[@]} -eq 0 ]; then
    echo "All datasets for $ADJ_FOLDER already processed!"
    exit 0
fi

# --- Run classifier on remaining CSVs ---
for CSV_FILE in "${TO_PROCESS[@]}"; do
    echo "Processing $CSV_FILE..."
    
    python "$CLASSIFIER_SCRIPT" "$CSV_FILE" "$OUTPUT_DIR"
    STATUS=$?

    if [ $STATUS -eq 0 ]; then
        echo "Classifier finished successfully for $CSV_FILE"
    else
        echo "ERROR: Classifier failed for $CSV_FILE" >&2
    fi
done

echo "✅ Finished all remaining datasets for $ADJ_FOLDER."
