#!/bin/bash

BOOK_CHAPTER_DIR="$HOME/confounded_analysis/scripts/evaluations/book_chapter"

# Cleanup old logs and temporary files to prevent space issues
echo "Cleaning up old logs and temporary files..."
find $HOME/.snakemake -type f -mtime +7 -delete 2>/dev/null || true

echo "Making directory"

OUTPUT_DIR="$HOME/confounded_analysis/grp_batch_effects/outputs/book_chapter"

mkdir -p "$OUTPUT_DIR"

echo "Changing permissions"
# For all subfolders, give all new files group ownership by default
find -L "$OUTPUT_DIR" -type d -exec chmod g+s {} + 

# --- NEW: Force Snakemake to find the CBC solver ---
# We ask pixi where cbc is, then export that path so Pulp sees it.
export PULP_CBC_PATH=$(pixi run which cbc)
echo "DEBUG: Manually setting CBC path to: $PULP_CBC_PATH"
# ---------------------------------------------------

SIMUL=2500

# Resources is tuned to the marylou cluster

echo "Starting Snakemake"
pixi run snakemake -s $BOOK_CHAPTER_DIR/Snakefile \
    --configfile $BOOK_CHAPTER_DIR/config.yaml \
    --max-jobs-per-second 10 \
    --max-status-checks-per-second 5 \
    --resources mem_mb=100000 runtime=4320 \
    --groups batch_simulation_group=20 \
    --groups batch_real_group=20 \
    --latency-wait 60 \
    --cores 28 \
    --rerun-incomplete --keep-going \
    --executor slurm --jobs $SIMUL \
    --envvars PATH CONDA_PREFIX

echo "Changing permissions"
# Give the group read and write access to all files and directories
find -L "$OUTPUT_DIR" -exec chmod g+wr {} +
# Give the group ownership over files and directories
chown -R :grp_batch_effects "$OUTPUT_DIR"