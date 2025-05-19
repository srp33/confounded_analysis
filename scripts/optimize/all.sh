#!/bin/bash

set -e

bash /scripts/optimize/run_scenarios.sh gse20194 unadjusted.csv batch treatment_response Sample "age,race,er_status,pr_status,bmn_grade,her2_status,histology,treatment_code"
bash /scripts/optimize/run_scenarios.sh gse24080 unadjusted.csv batch treatment_response Sample "age,race,er_status,pr_status,bmn_grade,her2_status,histology,treatment_code"
bash /scripts/optimize/run_scenarios.sh gse49711 unadjusted.csv Sex Class Sample_ID "INSS_Stage"

tmp_dir=/tmp/confounded

##################################################
# Summarizing optimization for gse20194
##################################################

out_dir=/outputs/optimizations/gse20194

mkdir -p ${out_dir}

combined_results_file=${out_dir}/combined_results.tsv
summarized_results_file=${out_dir}/summarized_results.tsv
graph_file=${out_dir}/summarized_results.pdf

python /scripts/optimize/cat_results.py "${tmp_dir}/results/gse20194/*" ${combined_results_file}
Rscript /scripts/optimize/summarize_results.R ${combined_results_file} unadjusted ${summarized_results_file}
Rscript /scripts/optimize/graph_results.R ${summarized_results_file} ${graph_file}

echo Summarized optimization results are saved in ${out_dir}

##################################################
# Summarizing optimization for gse49711
##################################################

out_dir=/outputs/optimizations/gse49711

mkdir -p ${out_dir}

combined_results_file=${out_dir}/combined_results.tsv
summarized_results_file=${out_dir}/summarized_results.tsv
graph_file=${out_dir}/summarized_results.pdf

python /scripts/optimize/cat_results.py "${tmp_dir}/results/gse49711/*" ${combined_results_file}
Rscript /scripts/optimize/summarize_results.R ${combined_results_file} unadjusted ${summarized_results_file}
Rscript /scripts/optimize/graph_results.R ${summarized_results_file} ${graph_file}

echo Summarized optimization results are saved in ${out_dir}

##################################################
# Combine the summarized results
##################################################

Rscript /scripts/optimize/combine_results.R /outputs/optimizations/gse20194/summarized_results.tsv /outputs/optimizations/gse49711/summarized_results.tsv /outputs/optimizations/combined_results.tsv /outputs/optimizations/combined_results.pdf
