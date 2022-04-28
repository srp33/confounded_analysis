#!/bin/bash

set -e

thisDir=$(dirname $0)

printf "\033[0;32mPreparing a simulated gene-expression dataset\033[0m\n"
python $thisDir/simulate_expression.py /data/simulated_expression/unadjusted.csv
