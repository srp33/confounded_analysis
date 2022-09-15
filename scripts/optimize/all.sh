#!/bin/bash

set -e

#bash /scripts/optimize/run_scenarios.sh simulated_expression unadjusted.csv Batch Class

#bash /scripts/optimize/run_scenarios.sh gse20194 unadjusted_small.csv treatment_response batch
bash /scripts/optimize/run_scenarios.sh gse20194 unadjusted.csv treatment_response batch

tmp_dir=/tmp/confounded

results_file=${tmp_dir}/combined_results.tsv

#python /scripts/optimize/cat_results.py "${tmp_dir}/results/gse20194/*" ${results_file}

#Rscript /scripts/optimize/summarize_results.R ${results_file} unadjusted
#Rscript /scripts/optimize/summarize_results.R ${results_file} unadjusted_small

#TODO: Try more hyperparameters on simulated data (see above).
#TODO: Do hyperparameter optimization with other datasets?
#        Need to tweak the way we compare against unadjusted metrics.
#TODO: Which parameter combinations are "good" for simulated data and real data?

#echo Results are saved in ${results_file}
