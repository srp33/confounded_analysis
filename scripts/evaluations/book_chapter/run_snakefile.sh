#!/bin/bash

BOOK_CHAPTER_DIR="$HOME/confounded_analysis/scripts/evaluations/book_chapter"
OUTPUT_DIR="$HOME/confounded_analysis/grp_batch_effects/outputs/book_chapter"

# Cleanup old logs and temporary files to prevent space issues
echo "Cleaning up old logs and temporary files..."
find $HOME/.snakemake -type f -mtime +7 -delete 2>/dev/null || true

# Setup output directory with group permissions
echo "Setting up output directory..."
mkdir -p "$OUTPUT_DIR"
find -L "$OUTPUT_DIR" -type d -exec chmod g+s {} + 


# The slurm account and slurm partition are essential for grouping
# Choosing the solver helps snakemake find it
# The resources and runtime are tuned to marylou, and help snakemake know what size of jobs it can schedule
echo "Starting Snakemake"
pixi run snakemake -s $BOOK_CHAPTER_DIR/Snakefile \
    --configfile $BOOK_CHAPTER_DIR/config.yaml \
    --scheduler-ilp-solver COIN_CMD \
    --executor slurm \
    --default-resources slurm_account=srp33 slurm_partition="(auto)" \
    --jobs 2500 \
    --group-components batch_real_group=20 batch_simulation_group=20 \
    --resources mem_mb=100000 runtime=4320 \
    --max-jobs-per-second 20 \
    --max-status-checks-per-second 10 \
    --latency-wait 120 \
    --rerun-incomplete \
    --keep-going