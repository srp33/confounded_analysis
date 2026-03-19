#!/usr/bin/env bash
#SBATCH --job-name=pathways
#SBATCH --output=logs/pathway_analysis/ttest_pathways_%A.log
#SBATCH --error=logs/pathway_analysis/ttest_pathways_%A.log
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G

set -euo pipefail

PATHWAY_SCRIPT="analyze_pathway.py"

# Directory where your gene list script wrote top genes CSV and .rnk files
GENE_LIST_DIR="/grphome/grp_batch_effects/outputs/metadata_features/target_pathways/gene_lists"

OUTDIR="/grphome/grp_batch_effects/outputs/metadata_features/target_pathways/ttest_analysis"
GMT_DIR="/grphome/grp_batch_effects/outputs/metadata_features/target_pathways/cancer_gmt"

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

module load python

mkdir -p "${OUTDIR}"

# Collect GMT files
mapfile -t GMT_FILES < <(find "${GMT_DIR}" -type f -name "*.gmt" ! -path "*pathway_analysis*" | sort)

if [[ ${#GMT_FILES[@]} -eq 0 ]]; then
    echo "ERROR: No GMT files found in ${GMT_DIR}"
    exit 1
fi

# Run pathway analysis — the Python script loops over each target's files
pixi run python "${PATHWAY_SCRIPT}" \
    --gene_lists_dir "${GENE_LIST_DIR}" \
    --gmt_files "${GMT_FILES[@]}" \
    --outdir "${OUTDIR}" \
    --top_n 100

echo "Done."