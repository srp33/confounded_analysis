#!/bin/bash

set -e

printf "\033[0;32mCalculating MMD\033[0m\n"

out_path="/outputs/metrics/mmd.csv"

rm -f ${out_path}

script_path="$(dirname $0)/mmd.py"

python "${script_path}" -i /data/simulated_expression -b Batch -o "$out_path"
python "${script_path}" -i /data/bladderbatch -b batch -o "$out_path"
python "${script_path}" -i /data/gse37199 -b plate -o "$out_path"
python "${script_path}" -i /data/tcga -b CancerType -o "$out_path"
python "${script_path}" -i /data/tcga_medium -b CancerType -o "$out_path"
python "${script_path}" -i /data/tcga_small -b CancerType -o "$out_path"
