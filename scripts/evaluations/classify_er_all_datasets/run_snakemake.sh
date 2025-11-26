#!/bin/bash
#SBATCH --job-name=run_snakemake
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=128
#SBATCH --mem=16G
#SBATCH --time=10:00:00
#SBATCH --output=logs/snakemake_%A_%a.out
#SBATCH --requeue  # allow job to be requeued if killed

pixi run snakemake -c 128




