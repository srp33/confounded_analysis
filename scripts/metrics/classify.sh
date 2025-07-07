#!/bin/bash

set -e

printf "\033[0;32mCalculating classification metrics\033[0m\n"

batch_out_path="/outputs/metrics/batch_classification.csv"
true_out_path="/outputs/metrics/true_classification.csv"

source /scripts/metrics/utils.sh

archive_file "${batch_out_path}"
archive_file "${true_out_path}"

# Remove lines pertaining to a particular adjuster and dataset
# sed -i '/wasserstein,gse20194,/d' "${batch_out_path}"
# sed -i '/wasserstein,gse20194,/d' "${true_out_path}"

script_path="$(dirname $0)/classify.py"


python "${script_path}" -i /data/gse20194 -o ${batch_out_path} -c meta_batch
python "${script_path}" -i /data/gse20194 -o ${true_out_path} -c meta_er_status
python "${script_path}" -i /data/gse20194 -o ${true_out_path} -c meta_her2_status
python "${script_path}" -i /data/gse20194 -o ${true_out_path} -c meta_pr_status

python "${script_path}" -i /data/gse24080 -o ${batch_out_path} -c meta_batch
python "${script_path}" -i /data/gse24080 -o ${true_out_path} -c meta_efs_outcome_label
python "${script_path}" -i /data/gse24080 -o ${true_out_path} -c meta_os_outcome_label


python "${script_path}" -i /data/gse49711 -o ${batch_out_path} -c meta_Class
# The classifier only does binary classification, so I split the classes multiple ways
python "${script_path}" -i /data/gse49711 -o ${true_out_path} -c meta_INSS_Stage_Split_3_4
