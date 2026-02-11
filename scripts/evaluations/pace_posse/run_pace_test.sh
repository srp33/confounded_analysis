#!/bin/bash

#SBATCH --time 6:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=4G
#SBATCH -J "pace_posse_workflow"
#SBATCH -o logs/snake_%A.log

# Pace/Posse independent workflow
OUTPUT_DIR="outputs"

# Cleanup old logs
echo "Cleaning up old logs..."
find $HOME/.snakemake -type f -mtime +7 -delete 2>/dev/null || true

# Setup output directory
echo "Setting up output directory..."
mkdir -p "$OUTPUT_DIR"
mkdir -p logs

echo "Starting Snakemake for PACE/POSSE workflow"
pixi run snakemake \
    --scheduler-ilp-solver COIN_CMD \
    --executor slurm \
    --default-resources slurm_account=srp33 slurm_partition="(auto)" \
    --jobs 100 \
    --resources mem_mb=50000 runtime=2160 \
    --max-jobs-per-second 10 \
    --max-status-checks-per-second 5 \
    --latency-wait 60 \
    --rerun-incomplete \
    2>&1 | tee logs/pace_posse_$(date +%Y%m%d_%H%M%S).log