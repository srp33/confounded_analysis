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

SIMUL=100

echo "Starting Snakemake"
snakemake -s $BOOK_CHAPTER_DIR/Snakefile \
    --configfile $BOOK_CHAPTER_DIR/config.yaml \
    --rerun-incomplete --printshellcmds --keep-going \
    --executor slurm --jobs $SIMUL \
    --envvars PATH CONDA_PREFIX

echo "Changing permissions"
# Give the group read and write access to all files and directories
find -L "$OUTPUT_DIR" -exec chmod g+wr {} +
# Give the group ownership over files and directories
chown -R :grp_batch_effects "$OUTPUT_DIR"