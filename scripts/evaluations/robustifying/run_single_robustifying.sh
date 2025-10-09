#!/bin/bash

# Wrapper script for running a single robustifying evaluation
# Parameters: N_sample_size max_batch_mean max_batch_var

if [ $# -ne 3 ]; then
    echo "Usage: $0 <N_sample_size> <max_batch_mean> <max_batch_var>"
    echo "Example: $0 20 3 4"
    exit 1
fi

N_SAMPLE_SIZE=$1
MAX_BATCH_MEAN=$2
MAX_BATCH_VAR=$3

echo "Starting robustifying evaluation with parameters:"
echo "  N_sample_size: $N_SAMPLE_SIZE"
echo "  max_batch_mean: $MAX_BATCH_MEAN"
echo "  max_batch_var: $MAX_BATCH_VAR"
echo "  Working directory: $(pwd)"
echo "  Date: $(date)"

# Change to the robustifying directory
cd /scripts/evaluations/robustifying

# Run the R script with the provided parameters
Rscript code/1_simpipe.R $N_SAMPLE_SIZE $MAX_BATCH_MEAN $MAX_BATCH_VAR

echo "Robustifying evaluation completed at $(date)"