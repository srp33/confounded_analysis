#!/bin/bash

#SBATCH --job-name=snakemake
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=128
#SBATCH --mem-per-cpu=2G
#SBATCH --time=04:00:00
#SBATCH --output=/grphome/grp_batch_effects/outputs/book_chapter/logs/snakemake_%j.out

BOOK_CHAPTER_DIR="/scripts/evaluations/book_chapter"

# Cleanup old logs and temporary files to prevent space issues
echo "Cleaning up old logs and temporary files..."
find $HOME/.snakemake -type f -mtime +7 -delete 2>/dev/null || true

mkdir -p $HOME/confounded_analysis/grp_batch_effects/outputs/book_chapter
chmod -R 777 $HOME/confounded_analysis/grp_batch_effects/outputs/book_chapter

SIMUL=$((SLURM_CPUS_PER_TASK-1))

bash $HOME/confounded_analysis/run_in_apptainer.sh snakemake -s $BOOK_CHAPTER_DIR/Snakefile --configfile $BOOK_CHAPTER_DIR/config.yaml -p --cores $SLURM_CPUS_PER_TASK --jobs $SIMUL --rerun-incomplete --printshellcmds --keep-going