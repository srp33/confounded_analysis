#!/bin/bash
#SBATCH --job-name=classify_unadjusted
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=10:00:00
#SBATCH --output=logs/classify_unadjusted_%A_%a.out
#SBATCH --requeue  # allow job to be requeued if killed

# --- Directories ---
CLASSIFIER_SCRIPT="/scripts/evaluations/classify_er_all_datasets/run_classifier.py"
DATA_DIR="/data/all_combined_subsets"
ADJUSTED_DIR="/data/adjusted_datasets"
OUTPUT_DIR="/outputs/classify_er_all"
UNADJUSTED_DIR="/data/adjusted_datasets/unadjusted"
OUT_SUBDIR="$OUTPUT_DIR/unadjusted"

mkdir -p "$OUT_SUBDIR"
mkdir -p "$UNADJUSTED_DIR"

# --- Generate Training/Testing sets ---
echo "Generating LOSO splits..."
python <<'PYCODE'
import pandas as pd
from pathlib import Path

data_dir = Path("/data/all_combined_subsets")
out_dir = Path("/data/adjusted_datasets/unadjusted")
out_dir.mkdir(parents=True, exist_ok=True)

for csv_path in sorted(data_dir.glob("subset_*studies.csv")):
    print(f"Creating train/test split for {csv_path}...")
    df = pd.read_csv(csv_path)
    if "meta_source" not in df.columns:
        raise ValueError(f"Missing 'meta_source' column in {csv_path}")

    studies = df["meta_source"].unique()
    base_name = csv_path.stem

    for study in studies:
        # The classifier script handles splitting internally
        out_path = out_dir / f"{base_name}_test_{study}.csv"
        df.to_csv(out_path, index=False)
        print(f"Created unadjusted combined file for test={study}: {out_path.name}")

PYCODE

# --- Get list of CSVs ---
shopt -s nullglob
CSV_FILES=("$UNADJUSTED_DIR"/*.csv)

if [ ${#CSV_FILES[@]} -eq 0 ]; then
    echo "No CSV files found in $UNADJUSTED_DIR"
    exit 1
fi 

# --- Run classifier on remaining CSVs ---
for CSV_FILE in "${CSV_FILES[@]}"; do
    echo "Processing $CSV_FILE..."
    
    python "$CLASSIFIER_SCRIPT" "$CSV_FILE" "$OUT_SUBDIR"
    STATUS=$?

    if [ $STATUS -eq 0 ]; then
        echo "Classifier finished successfully for $CSV_FILE"
    else
        echo "ERROR: Classifier failed for $CSV_FILE" >&2
    fi
done

echo "✅ Finished classifying all unadjusted datasets."
