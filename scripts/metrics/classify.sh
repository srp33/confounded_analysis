#!/bin/bash

set -e

printf "\033[0;32mCalculating classification metrics\033[0m\n"

batch_out_path="/outputs/metrics/batch_classification.csv"
true_out_path="/outputs/metrics/true_classification.csv"

rm -f ${batch_out_path} ${true_out_path}

script_path="$(dirname $0)/classify.py"

# python "${script_path}" -i /data/simulated_expression -o ${batch_out_path} -c Batch 
# python "${script_path}" -i /data/simulated_expression -o ${true_out_path} -c Class


# It doesn't make sense to measure true-class accuracy for bladderbatch because it is confounded with batch
python "${script_path}" -i /data/bladderbatch -o ${batch_out_path} -c batch

# python "${script_path}" -i /data/gse37199 -o ${batch_out_path} -c plate
# python "${script_path}" -i /data/gse37199 -o ${true_out_path} -c Stage

for dataset in tcga tcga_medium tcga_small
do
  python "${script_path}" -i /data/$dataset -o ${batch_out_path} -c CancerType
  python "${script_path}" -i /data/$dataset -o ${true_out_path} -c TP53_Mutated
done
