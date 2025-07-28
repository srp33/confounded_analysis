#!/bin/bash

set -e

printf "\033[0;32mCalculating classification metrics\033[0m\n"

batch_out_path="/outputs/metrics/batch_classification.csv"
true_out_path="/outputs/metrics/true_classification.csv"

source /scripts/metrics/utils.sh

archive_file "${batch_out_path}"
archive_file "${true_out_path}"

# Remove lines pertaining to a particular adjuster and dataset
# sed -i '/npn,gse20194,/d' "${batch_out_path}"
# sed -i '/npn,gse20194,/d' "${true_out_path}"

script_path="$(dirname $0)/classify.py"



# Predict, globally
batch_out_path_global="/outputs/metrics/global_batch_classification.csv"
true_out_path_global="/outputs/metrics/global_true_classification.csv"

python "${script_path}" -i /data/gse20194 -o ${batch_out_path_global} -p meta_batch
python "${script_path}" -i /data/gse20194 -o ${true_out_path_global} -p meta_er_status

python "${script_path}" -i /data/gse20194 -o ${batch_out_path_global} -p meta_batch
python "${script_path}" -i /data/gse20194 -o ${true_out_path_global} -p meta_her2_status

python "${script_path}" -i /data/gse20194 -o ${batch_out_path_global} -p meta_batch
python "${script_path}" -i /data/gse20194 -o ${true_out_path_global} -p meta_pr_status

# 24080
python "${script_path}" -i /data/gse24080 -o ${batch_out_path_global} -p meta_batch
python "${script_path}" -i /data/gse24080 -o ${true_out_path_global} -p meta_cytogenetic_abnormality


# 49711
# The classifier only does binary classification, so I split the classes multiple ways
python "${script_path}" -i /data/gse49711 -o ${batch_out_path_global} -p meta_Class
python "${script_path}" -i /data/gse49711 -o ${true_out_path_global} -p meta_INSS_Stage_Split_3_4



# Predict, conditional on the other column

# 20194
python "${script_path}" -i /data/gse20194 -o ${batch_out_path} -p meta_batch -c meta_er_status
python "${script_path}" -i /data/gse20194 -o ${true_out_path} -p meta_er_status -c meta_batch

python "${script_path}" -i /data/gse20194 -o ${batch_out_path} -p meta_batch -c meta_her2_status
python "${script_path}" -i /data/gse20194 -o ${true_out_path} -p meta_her2_status -c meta_batch

python "${script_path}" -i /data/gse20194 -o ${true_out_path} -p meta_pr_status -c meta_batch
python "${script_path}" -i /data/gse20194 -o ${batch_out_path} -p meta_batch -c meta_pr_status

# 24080
python "${script_path}" -i /data/gse24080 -o ${true_out_path} -p meta_cytogenetic_abnormality -c meta_batch
python "${script_path}" -i /data/gse24080 -o ${batch_out_path} -p meta_batch -c meta_cytogenetic_abnormality


# 49711
# The classifier only does binary classification, so I split the classes multiple ways
python "${script_path}" -i /data/gse49711 -o ${true_out_path} -p meta_INSS_Stage_Split_3_4 -c meta_Class
python "${script_path}" -i /data/gse49711 -o ${batch_out_path} -p meta_Class -c meta_INSS_Stage_Split_3_4

