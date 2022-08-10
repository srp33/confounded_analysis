#!/bin/bash

set -e

bash /scripts/optimize/run_scenarios.sh simulated_expression Batch

#rm -f ${results_file}
#python /scripts/optimize/cat_results.py "${tmp_dir}/results/*" ${results_file}

#Rscript scripts/optimize/summarize_results.R ${results_file}

#TODO: Try more hyperparameters on simulated data (see above).
#TODO: Do hyperparameter optimization with other datasets?
#        Need to tweak the way we compare against unadjusted metrics.
#TODO: Which parameter combinations are "good" for simulated data and real data?

#echo Results are saved in ${results_file}
