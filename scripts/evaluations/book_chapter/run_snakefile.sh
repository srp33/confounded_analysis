#!/bin/bash

#SBATCH --job-name=snakemake
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH --mem-per-cpu=2G
#SBATCH --time=02:00:00
#SBATCH --output=snakemake_%j.out

BOOK_CHAPTER_DIR="/scripts/evaluations/book_chapter"

mkdir -p $HOME/confounded_analysis/grp_batch_effects/outputs/book_chapter
chmod -R 777 $HOME/confounded_analysis/grp_batch_effects/outputs/book_chapter

SIMUL=$((SLURM_CPUS_PER_TASK-2))

bash $HOME/confounded_analysis/run_in_apptainer.sh snakemake -s $BOOK_CHAPTER_DIR/Snakefile --configfile $BOOK_CHAPTER_DIR/config.yaml -p --cores $SLURM_CPUS_PER_TASK --jobs $SIMUL --rerun-incomplete --printshellcmds --keep-going