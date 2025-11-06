#!/bin/bash

#SBATCH --job-name=snakemake
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=256G
#SBATCH --time=06:00:00
#SBATCH --output=snakemake_%j.out

BOOK_CHAPTER_DIR="/scripts/evaluations/book_chapter"

mkdir -p $HOME/confounded_analysis/grp_batch_effects/outputs/book_chapter
chmod -R 777 $HOME/confounded_analysis/grp_batch_effects/outputs/book_chapter

# 128 cores, up to 112 simultaneous R processes
bash $HOME/confounded_analysis/run_in_apptainer.sh snakemake -s $BOOK_CHAPTER_DIR/Snakefile --configfile $BOOK_CHAPTER_DIR/config.yaml -p --cores 128 --jobs 112 --rerun-incomplete --printshellcmds --keep-going