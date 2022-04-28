#!/bin/bash

set -e

printf "\033[0;32mPreparing a simulated gene-expression dataset\033[0m\n"
python simulate_expression.py /data/simulated_expression/unadjusted.csv
