#!/bin/bash

# Wrapper script to run 1_simpipe.R with parameters
# Arguments: N_sample_size, max_batch_mean, max_batch_var

if [ $# -ne 3 ]; then
    echo "Usage: $0 <N_sample_size> <max_batch_mean> <max_batch_var>"
    exit 1
fi

N_SAMPLE_SIZE=$1
MAX_BATCH_MEAN=$2
MAX_BATCH_VAR=$3

echo "Running robustifying analysis with parameters:"
echo "  N_sample_size: $N_SAMPLE_SIZE"
echo "  max_batch_mean: $MAX_BATCH_MEAN"
echo "  max_batch_var: $MAX_BATCH_VAR"

# Change to the robustifying directory and run the R script
cd /scripts/evaluations/robustifying
Rscript code/1_simpipe.R $N_SAMPLE_SIZE $MAX_BATCH_MEAN $MAX_BATCH_VAR