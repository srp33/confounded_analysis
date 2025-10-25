#!/bin/bash
#SBATCH --job-name=build_combat_container
#SBATCH --qos=login
#SBATCH --partition=login
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=build_combat_apptainer.out

module load apptainer

echo "Checking internet connection..."
if ! wget -q --spider http://google.com; then
    echo "ERROR: No internet connection. Apptainer build requires internet access." >&2
    exit 1
fi

echo "Starting Apptainer build at $(date)"
apptainer build \
    ~/confounded_analysis/apptainer/combatseq_image.sif \
    ~/confounded_analysis/apptainer/apptainer_combatseq.def
echo "Build finished at $(date)"