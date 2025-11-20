#!/bin/bash

# sbatch --mem=4G --time=23:59:59 -c 2 -o snake.out run_snakefile.sh

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

SIMUL=2500

echo "Starting Snakemake"
pixi run snakemake -s $BOOK_CHAPTER_DIR/Snakefile \
    --configfile $BOOK_CHAPTER_DIR/config.yaml \
    --max-jobs-per-second 5 \
    --max-status-checks-per-second 10 \
    --groups batch_simulation_group=20 \
    --groups batch_real_group=20 \
    --latency-wait 60 \
    --rerun-incomplete --keep-going \
    --executor slurm --jobs $SIMUL \
    --envvars PATH CONDA_PREFIX

echo "Changing permissions"
# Give the group read and write access to all files and directories
find -L "$OUTPUT_DIR" -exec chmod g+wr {} +
# Give the group ownership over files and directories
chown -R :grp_batch_effects "$OUTPUT_DIR"