#!/bin/bash
#SBATCH --job-name=tstats_batch
#SBATCH --output=logs/ttest/tstats_%j.log
#SBATCH --error=logs/ttest/tstats_%j.log
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G

ADJUST_DIR="/grphome/grp_batch_effects/outputs/metadata_features/labeled_adjusted"
OUTDIR="/grphome/grp_batch_effects/outputs/metadata_features/ttest"

mkdir -p "$OUTDIR"
TARGETS=(
    meta_er_status
    meta_her2_status
    meta_sex
    meta_chemotherapy
    meta_age_at_diagnosis_combined_lt50
    meta_age_at_diagnosis_combined_50_69
    meta_age_at_diagnosis_combined_ge70
    meta_menopause_status
    meta_histological_type
    )

# Collect all adjusted CSVs
mapfile -t CSV_FILES < <(find "${ADJUST_DIR}" -name "*.csv" | sort)

python t_test.py \
    "${CSV_FILES[@]}" \
    --meta-cols "${TARGETS[@]}" \
    --outdir "${OUTDIR}"

echo "Finished at $(date)"