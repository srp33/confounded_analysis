#!/bin/bash

set -e

###bash /scripts/optimize/run_scenarios.sh simulated_expression unadjusted.csv Batch Class

#bash /scripts/optimize/run_scenarios.sh gse20194 unadjusted.csv batch treatment_response Sample "age,race,er_status,pr_status,bmn_grade,her2_status,histology,treatment_code"

tmp_dir=/tmp/confounded
out_dir=/outputs/optimizations/gse20194

mkdir -p ${out_dir}

combined_results_file=${out_dir}/combined_results.tsv
summarized_results_file=${out_dir}/summarized_results.tsv

python /scripts/optimize/cat_results.py "${tmp_dir}/results/gse20194/*" ${combined_results_file}
Rscript /scripts/optimize/summarize_results.R ${combined_results_file} unadjusted ${summarized_results_file}
Rscript /scripts/optimize/graph_results.R ${summarized_results_file} 

echo Summarized optimization results are saved in ${out_dir}
