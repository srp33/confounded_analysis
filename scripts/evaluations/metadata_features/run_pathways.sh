#!/usr/bin/env bash
#SBATCH --job-name=pathways
#SBATCH --output=logs/pathway_analysis/pathways_%A_%a.log
#SBATCH --error=logs/pathway_analysis/pathways_%A_%a.log
#SBATCH --time=16:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G

set -euo pipefail

PATHWAY_SCRIPT="analyze_pathway.py"
PERM_DIR="/grphome/grp_batch_effects/outputs/metadata_features/permutation_importance"
OUTDIR="/grphome/grp_batch_effects/outputs/metadata_features/pathway_analysis"
META_DIR="/grphome/grp_batch_effects/outputs/metadata_features"
SELECTED_GENES_CSV="/grphome/grp_batch_effects/outputs/metadata_features/plots/selected_genes.csv"

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

module load python

mkdir -p "${OUTDIR}"

# -------------------------------------------------
# Collect GMT Files
# -------------------------------------------------

mapfile -t GMT_FILES < <(find "${META_DIR}" -type f -name "*.gmt" ! -path "*pathway_analysis*" | sort)

if [[ ${#GMT_FILES[@]} -eq 0 ]]; then
    echo "ERROR: No GMT files found in ${META_DIR}"
    exit 1
fi

pixi run python "${PATHWAY_SCRIPT}" \
  --selected_genes_csv "${SELECTED_GENES_CSV}" \
  --gmt_files "${GMT_FILES[@]}" \
  --outdir "${OUTDIR}"

echo "Done."
