#!/bin/bash
#SBATCH --job-name=classify_all
#SBATCH --array=0 # One job per adjuster folder
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=logs/classify_%A_%a.out

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
ANALYSIS_DIR=~/confounded_analysis
CLASSIFIER_SCRIPT=/scripts/evaluations/classify_er_all_datasets/run_classifier.sh
UNADJUSTED_SCRIPT=/scripts/evaluations/classify_er_all_datasets/classify_unadjusted.sh

echo "Running classification with $CLASSIFIER_SCRIPT"
echo "Running classification on unadjusted files with $UNADJUSTED_SCRIPT"
echo "Job ID: $SLURM_JOB_ID"

cd $ANALYSIS_DIR
#bash $ANALYSIS_DIR/run_in_apptainer.sh $CLASSIFIER_SCRIPT "$SLURM_ARRAY_TASK_ID"

if [ "$SLURM_ARRAY_TASK_ID" -eq 0 ]; then
    echo "Running classification on unadjusted files with $UNADJUSTED_SCRIPT"
    bash $ANALYSIS_DIR/run_in_apptainer.sh $UNADJUSTED_SCRIPT
else 
    echo "Skipping unadjusted classification for task $SLURM_ARRAY_TASK_ID"
fi