#!/bin/bash

set -e

Rscript scripts/figures/er_classification_plots_single.R
bash scripts/figures/classification_figures.sh
bash scripts/figures/tsne.sh
bash scripts/figures/plot_reduced.sh
bash scripts/figures/mse_mmd_classification.sh
