#!/bin/bash

set -e

printf "\033[0;32mGenerating Classification figures\033[0m\n"

# Partial
Rscript --vanilla  scripts/figures/classification_figures.R -i "/outputs/metrics" -f "/outputs/figures/classification/" -s "/../data/" -b "batch_classification.csv" -t "true_classification.csv"

# Global
Rscript --vanilla  scripts/figures/classification_figures.R -i "/outputs/metrics" -f "/outputs/figures/classification_global/" -s "/../data/" -b "global_batch_classification.csv" -t "global_true_classification.csv"
