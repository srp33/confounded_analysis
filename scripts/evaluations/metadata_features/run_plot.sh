#!/usr/bin/env bash
#SBATCH --job-name=perm-heatmaps
#SBATCH --output=logs/plot/heatmaps_%j.log
#SBATCH --error=logs/plot/heatmaps_%j.log
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
set -euo pipefail

# -------------------------
# User-configurable paths
# -------------------------
PLOT_SCRIPT="plot_heatmap.py"
PERM_DIR="/grphome/grp_batch_effects/outputs/metadata_features/permutation_importance"
OUTDIR="/grphome/grp_batch_effects/outputs/metadata_features/plots"
TRAIN="gse62944_tumor"
TEST="metabric"
ALIGNED_CSV="/grphome/grp_batch_effects/outputs/metadata_features/subset/aligned_subset.csv"

# Optional parameters
THRESHOLD=0.004

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

module load python

# -------------------------
# Collect CSVs (one per adjuster)
# -------------------------
CSV_FILES=$(find "${PERM_DIR}" -type f -name "2_*_permutation_importance.csv" | sort)

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
  --threshold "${THRESHOLD}" \
  --train "${TRAIN}" \
  --test "${TEST}" \
  --aligned_csv "${ALIGNED_CSV}"

echo
echo "All files written to: ${OUTDIR}"