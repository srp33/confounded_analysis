#!/bin/bash

set -e

printf "\033[0;32mGenerating Classification figures\033[0m\n"

# Global
# Rscript --vanilla  scripts/figures/classification_figures.R -i "/outputs/metrics" -f "/outputs/figures/classification_global/" -s "/../data/gold/" -b "global_batch_classification.csv" -t "global_true_classification.csv"

Rscript --vanilla  scripts/figures/combo_classification_figures.R -i "/outputs/metrics" -f "/outputs/figures/combined_global/" -s "/../data/combined_data/" -b "combined_global_batch_classification.csv" -t "combined_global_true_classification.csv"


# Partial
# Rscript --vanilla  scripts/figures/classification_figures.R -i "/outputs/metrics" -f "/outputs/figures/classification/" -s "/../data/gold/" -b "batch_classification.csv" -t "true_classification.csv"

Rscript --vanilla  scripts/figures/combo_classification_figures.R -i "/outputs/metrics" -f "/outputs/figures/combined/" -s "/../data/combined_data/" -b "combined_batch_classification.csv" -t "combined_true_classification.csv"
