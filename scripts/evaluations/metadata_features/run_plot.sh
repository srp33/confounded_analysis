#!/usr/bin/env bash
#SBATCH --job-name=perm-heatmaps
#SBATCH --output=logs/plot/heatmaps_%j.log
#SBATCH --error=logs/plot/heatmaps_%j.log
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
set -euo pipefail

# -------------------------
# User-configurable paths
# -------------------------
PLOT_SCRIPT="plot_heatmap.py"
PERM_DIR="/grphome/grp_batch_effects/outputs/metadata_features/permutation_importance"
OUTDIR="/grphome/grp_batch_effects/outputs/metadata_features/plots"

# Optional parameters
THRESHOLD=0.003

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

module load python

# -------------------------
# Collect CSVs (one per adjuster)
# -------------------------
CSV_FILES=$(find "${PERM_DIR}" -type f -name "*_permutation_importance.csv" | sort)

if [[ -z "${CSV_FILES}" ]]; then
  echo "ERROR: No permutation importance CSVs found in ${PERM_DIR}"
  exit 1
fi

echo "Found permutation importance CSVs:"
echo "${CSV_FILES}"
echo

# -------------------------
# Create output directory
# -------------------------
mkdir -p "${OUTDIR}"

# -------------------------
# Run heatmap plotting
# -------------------------
pixi run python "${PLOT_SCRIPT}" \
  --csvs ${CSV_FILES} \
  --outdir "${OUTDIR}" \
  --threshold "${THRESHOLD}"

echo
echo "All heatmaps written to: ${OUTDIR}"