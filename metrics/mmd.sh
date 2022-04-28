#!/bin/bash

set -e

printf "\033[0;32mCalculating MMD\033[0m\n"

out_path="/output/metrics/mmd.csv"

python mmd.py -i /data/simulated_expression -b Batch -o $out_path
#python mmd.py -i /data/mnist -b Batch -o $out_path
#python mmd.py -i /data/bladderbatch -b batch -o $out_path
#python mmd.py -i /data/gse37199 -b plate -o $out_path
#python mmd.py -i /data/tcga -b CancerType -o $out_path
#python mmd.py -i /data/tcga_medium -b CancerType -o $out_path
#python mmd.py -i /data/tcga_small -b CancerType -o $out_path
