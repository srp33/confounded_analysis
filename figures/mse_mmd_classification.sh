#!/bin/bash

set -e

printf "\033[0;32mGenerating MSE, MMD, and classification figures and tables\033[0m\n"

Rscript --vanilla mse_mmd_classification.R
