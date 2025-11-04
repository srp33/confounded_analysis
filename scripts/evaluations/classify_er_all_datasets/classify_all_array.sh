#!/bin/bash
#SBATCH --job-name=classify_all
#SBATCH --array=1-4
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=logs/adjuster_%A_%a.out

module load apptainer

ADJUSTERS=("min_mean" "log_combat" "mnn" "gmm")
ADJ=${ADJUSTERS[$SLURM_ARRAY_TASK_ID-1]}

echo "Running adjuster: $ADJ"

apptainer exec ~/confounded_analysis/apptainer/remove-batch-effects.sif Rscript scripts/evaluations/classify_er_all_datasets/run_scaling_experiment.R --adjuster $ADJ