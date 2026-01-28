#!/bin/bash
#SBATCH --job-name=adjust_data
#SBATCH --array=0-4 # 5 adjusters
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=logs/adjust/adjuster_%A_%a.out

# Set up fixed inputs
SCALING_EXPERIMENT_FILE="../classify_er_all_datasets/run_scaling_experiment.R"
SUBSET_PATH="/grphome/grp_batch_effects/outputs/metadata_features/subset/combined.csv"
K=2
TEST_SOURCE="metabric"
OUTPUT_DIR="/grphome/grp_batch_effects/outputs/metadata_features/adjusted"
ADJUST_SCRIPT="../../adjust/adjust.R"
METADATA_FILE="/grphome/grp_batch_effects/data/geo_metadata.csv"

# Adjusters to run
ADJUSTERS=("mnn" "min_mean" "log_combat" "gmm" "log_transformed")

# Pick adjuster for this task
ADJ="${ADJUSTERS[$SLURM_ARRAY_TASK_ID]}"

echo "SLURM job: ${SLURM_JOB_ID}"
echo "Array task: ${SLURM_ARRAY_TASK_ID}"
echo "Running adjuster: ${ADJ}"

# Run one adjuster
pixi run Rscript ${SCALING_EXPERIMENT_FILE} \
    --adjuster "${ADJ}" \
    --subset-path "${SUBSET_PATH}" \
    --k "${K}" \
    --test "${TEST_SOURCE}" \
    --output-dir "${OUTPUT_DIR}" \
    --adjust-script "${ADJUST_SCRIPT}" \
    --metadata-file "${METADATA_FILE}"

echo "Finished adjuster: ${ADJ}."