#!/bin/bash

set -e

printf "\033[0;32mGenerating Classification figures\033[0m\n"

Rscript --vanilla  scripts/figures/classification_figures.R