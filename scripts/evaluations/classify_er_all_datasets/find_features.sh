#!/bin/bash
#SBATCH --job-name=features
#SBATCH --output=log/features_%j.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=10:00:00

# -----------------------------
# CONFIGURATION
# -----------------------------

# Directory containing CSV results
FEATURE_SCRIPT="classify_feature_importance.py"
DATA_DIR="/grphome/grp_batch_effects/data/adjusted_data"
OUT_DIR="/grphome/grp_batch_effects/data/feature_importance"

# ----------------------------
# Loop through adjusted data
# ----------------------------
mapfile -t CSV_FILES < <(find "$DATA_DIR" -type f -name "*.csv")

if [ ${#CSV_FILES[@]} -eq 0 ]; then
    echo "No CSV files found under $DATA_DIR"
    exit 1
fi 

echo "Found ${#CSV_FILES[@]} CSV files"

# -----------------------------
# Classify feature importance
# -----------------------------
for CSV_FILE in "${CSV_FILES[@]}"; do
    echo "Processing $CSV_FILE..."

    pixi run python "$FEATURE_SCRIPT" --csv "$CSV_FILE" --outdir "$OUT_DIR"
    STATUS=$?

    if [ "$STATUS" -eq 0 ]; then
        echo "Classifier finished successfully for $CSV_FILE"
    else 
        echo "ERROR: Classifier failed for $CSV_FILE" >&2
    fi 
done

echo "Finished all remaining datasets for $DATA_DIR"