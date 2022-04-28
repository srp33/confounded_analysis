#!/bin/bash

set -e

printf "\033[0;32mGenerating loss figure\033[0m\n"

Rscript --vanilla loss.R
