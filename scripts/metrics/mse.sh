#!/bin/bash

set -e

printf "\033[0;32mCalculating MSE\033[0m\n"

out_path="/outputs/metrics/mse.csv"

rm -f ${out_path}

script_path="$(dirname $0)/mse.py"

python "${script_path}" -i /data/simulated_expression -o "$out_path"
python "${script_path}" -i /data/bladderbatch -o "$out_path"
python "${script_path}" -i /data/gse37199 -o "$out_path"
python "${script_path}" -i /data/tcga -o "$out_path"
python "${script_path}" -i /data/tcga_medium -o "$out_path"
python "${script_path}" -i /data/tcga_small -o "$out_path"
