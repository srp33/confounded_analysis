#!/bin/bash

set -e

printf "\033[0;32mCalculating MSE\033[0m\n"

out_path="/output/metrics/mse.csv"

python mse.py -i /data/simulated_expression -o $out_path
#python mse.py -i /data/mnist -o $out_path
#python mse.py -i /data/bladderbatch -o $out_path
#python mse.py -i /data/gse37199 -o $out_path
#python mse.py -i /data/tcga -o $out_path
#python mse.py -i /data/tcga_medium -o $out_path
#python mse.py -i /data/tcga_small -o $out_path
