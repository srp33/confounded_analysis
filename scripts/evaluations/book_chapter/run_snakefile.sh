#!/bin/bash

#SBATCH --job-name=snakemake
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=128
#SBATCH --mem-per-cpu=2G
#SBATCH --time=04:00:00
#SBATCH --output=/grphome/grp_batch_effects/outputs/book_chapter/logs/snakemake_%j.out

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

SLURM_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK:-16}
SIMUL=$((SLURM_CPUS_PER_TASK-1))

echo "SLURM_CPUS_PER_TASK $SLURM_CPUS_PER_TASK"
echo  "SIMUL $SIMUL"

echo "Starting Snakemake"
snakemake -s $BOOK_CHAPTER_DIR/Snakefile --configfile $BOOK_CHAPTER_DIR/config.yaml -p --cores $SLURM_CPUS_PER_TASK --jobs $SIMUL --rerun-incomplete --printshellcmds --keep-going

echo "Changing permissions"
# Give the group read and write access to all files and directories
find -L "$OUTPUT_DIR" -exec chmod g+wr {} +
# Give the group ownership over files and directories
chown -R :grp_batch_effects "$OUTPUT_DIR"