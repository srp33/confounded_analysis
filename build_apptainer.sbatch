#!/bin/bash
#SBATCH --job-name=build_container
#SBATCH --qos=login
#SBATCH --partition=login
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=build_apptainer.out
#SBATCH --error=build_apptainer.err

module load apptainer

echo "Starting Apptainer build at $(date)"
apptainer build \
    ~/confounded_analysis/apptainer/remove-batch-effects.sif \
    ~/confounded_analysis/apptainer/apptainer.def
echo "Build finished at $(date)"