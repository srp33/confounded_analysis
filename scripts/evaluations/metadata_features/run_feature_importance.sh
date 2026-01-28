#!/bin/bash
#SBATCH --job-name=feature_importance
#SBATCH --array=0-4           # one job per adjuster
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=logs/feature_importance/fi_%A_%a.out

set -euo pipefail

# -------------------------
# Paths
# -------------------------
FEATURE_SCRIPT="../classify_er_all_datasets/classify_feature_importance.py"
CLASSIFY_DIR="/grphome/grp_batch_effects/outputs/metadata_features/classify"
OUTPUT_DIR="/grphome/grp_batch_effects/outputs/metadata_features/feature_importance"

# -------------------------
# Adjusters
# -------------------------
ADJUSTERS=("mnn" "min_mean" "log_combat" "gmm" "log_transformed")

# -------------------------
# Select adjuster for this job
# -------------------------
ADJ="${ADJUSTERS[$SLURM_ARRAY_TASK_ID]}"
echo "======================================"
echo "SLURM job: $SLURM_JOB_ID"
echo "Task ID: $SLURM_ARRAY_TASK_ID"
echo "Processing adjuster: $ADJ"
echo "======================================"

# -------------------------
# Find all classify CSVs for this adjuster
# -------------------------
mapfile -t CSV_FILES < <(find "${CLASSIFY_DIR}/${ADJ}" -name "*.csv" | sort)

if [ ${#CSV_FILES[@]} -eq 0 ]; then
    echo "No CSV files found for adjuster ${ADJ}"
    exit 1
fi

# -------------------------
# Aggregate CSVs
# -------------------------
AGGREGATED_CSV="${OUTPUT_DIR}/${ADJ}_all_chunks.csv"
mkdir -p "${OUTPUT_DIR}"

# Use header from the first file
head -n 1 "${CSV_FILES[0]}" > "${AGGREGATED_CSV}"

# Append data from all CSVs
for f in "${CSV_FILES[@]}"; do
    tail -n +2 "$f" >> "${AGGREGATED_CSV}"
done

echo "Aggregated ${#CSV_FILES[@]} files into ${AGGREGATED_CSV}"

# -------------------------
# Run feature importance
# -------------------------
pixi run python "${FEATURE_SCRIPT}" \
    --csv "${AGGREGATED_CSV}" \
    --outdir "${OUTPUT_DIR}" \
    --metric "roc_auc" \
    --n_jobs 4

echo "Feature importance completed for ${ADJ}"
