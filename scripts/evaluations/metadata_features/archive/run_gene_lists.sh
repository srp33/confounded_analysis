#!/usr/bin/env bash
#SBATCH --job-name=gene_lists
#SBATCH --output=logs/gene_list/gene_list_%j.log
#SBATCH --error=logs/gene_list/gene_list_%j.log
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
set -euo pipefail

# -------------------------
# User-configurable paths
# -------------------------
GENE_LIST_SCRIPT="generate_gene_lists.py"
PLOT_SCRIPT="plot.py"
#PERM_DIR="/grphome/grp_batch_effects/outputs/metadata_features/permutation_importance"
TTEST_DIR="/grphome/grp_batch_effects/outputs/metadata_features/ttest"
OUTDIR="/grphome/grp_batch_effects/outputs/metadata_features/target_pathways"
GENE_LIST_DIR="${OUTDIR}/gene_lists"

TRAIN="gse62944_tumor"
TEST="metabric"
ALIGNED_CSV="/grphome/grp_batch_effects/outputs/metadata_features/subset/aligned_subset.csv"

# Optional parameters
THRESHOLD=0.005
TOP_N=100

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

module load python

# -------------------------
# Collect CSVs (one per adjuster)
# -------------------------
CSV_FILES=$(find "${TTEST_DIR}" -type f -name "*.csv" | sort)

if [[ -z "${CSV_FILES}" ]]; then
  echo "ERROR: No ttest CSVs found in ${TTEST_DIR}"
  exit 1
fi

echo "Found ttest CSVs:"
echo "${CSV_FILES}"
echo

# -------------------------
# Create output directory
# -------------------------
mkdir -p "${OUTDIR}"
mkdir -p "${GENE_LIST_DIR}"

# -------------------------
# Run gene list
# -------------------------
pixi run python "${GENE_LIST_SCRIPT}" \
  --csvs ${CSV_FILES} \
  --outdir ${GENE_LIST_DIR} \
  --threshold ${THRESHOLD}

# # -------------------------
# # Collect ranked list files
# # -------------------------
# RANKED_LISTS=$(find "${GENE_LIST_DIR}" -type f -name "*_ranked.csv" | sort)

# if [[ -z "${RANKED_LISTS}" ]]; then 
#     echo "ERROR: No ranked gene lists found."
#     exit 1
# fi 

# # -------------------------
# # Run plotting 
# # -------------------------
# pixi run python "${PLOT_SCRIPT}" \
#     --expression_csv "${ALIGNED_CSV}" \
#     --ranked_lists ${RANKED_LISTS} \
#     --outdir "${OUTDIR}" \
#     --top_n "${TOP_N}"

# echo
# echo "All files written to: ${OUTDIR}"