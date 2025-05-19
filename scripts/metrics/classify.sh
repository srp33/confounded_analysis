#!/bin/bash

set -e

printf "\033[0;32mCalculating classification metrics\033[0m\n"

batch_out_path="/outputs/metrics/batch_classification.csv"
true_out_path="/outputs/metrics/true_classification.csv"

rm -f ${batch_out_path} ${true_out_path}

script_path="$(dirname $0)/classify.py"

# python "${script_path}" -i /data/simulated_expression -o ${batch_out_path} -c Batch 
# python "${script_path}" -i /data/simulated_expression -o ${true_out_path} -c Class


python "${script_path}" -i /data/gse20194 -o ${batch_out_path} -c batch
python "${script_path}" -i /data/gse20194 -o ${true_out_path} -c er_status
python "${script_path}" -i /data/gse20194 -o ${true_out_path} -c her2_status
python "${script_path}" -i /data/gse20194 -o ${true_out_path} -c pr_status

python "${script_path}" -i /data/gse24080 -o ${batch_out_path} -c batch
python "${script_path}" -i /data/gse24080 -o ${true_out_path} -c efs_outcome_label
python "${script_path}" -i /data/gse24080 -o ${true_out_path} -c os_outcome_label

python "${script_path}" -i /data/gse49711 -o ${batch_out_path} -c Class
python "${script_path}" -i /data/gse49711 -o ${true_out_path} -c INSS_Stage