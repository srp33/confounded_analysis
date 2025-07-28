#!/bin/bash

set -e

bash scripts/figures/plot_reduced.sh
bash scripts/figures/classification_figures.sh
# bash scripts/figures/loss.sh
# bash scripts/figures/tsne.sh
bash scripts/figures/mse_mmd_classification.sh
