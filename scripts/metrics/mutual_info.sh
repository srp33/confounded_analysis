#!/bin/bash

set -e

printf "\033[0;32mCalculating Mutual info\033[0m\n"

out_path="/outputs/metrics/mutual_info.csv"
confusion_path="/outputs/metrics/confusion_matrix.md"

source /scripts/metrics/utils.sh

# Save previous file to an archive
archive_file "${out_path}"
archive_file "${confusion_path}"

script_path="$(dirname $0)/mutual_info.py"

python "${script_path}" -i /data/gse20194 -o ${out_path} -c meta_er_status -b meta_batch -m ${confusion_path}
python "${script_path}" -i /data/gse20194 -o ${out_path} -c meta_her2_status -b meta_batch -m ${confusion_path}
python "${script_path}" -i /data/gse20194 -o ${out_path} -c meta_pr_status -b meta_batch -m ${confusion_path}

python "${script_path}" -i /data/gse24080 -o ${out_path} -c meta_efs_outcome_label -b meta_batch -m ${confusion_path}
python "${script_path}" -i /data/gse24080 -o ${out_path} -c meta_os_outcome_label -b meta_batch -m ${confusion_path}

# The classifier only does binary classification, so I split the classes multiple ways
python "${script_path}" -i /data/gse49711 -o ${out_path} -c meta_INSS_Stage_Split_1_2 -b meta_Class -m ${confusion_path}
python "${script_path}" -i /data/gse49711 -o ${out_path} -c meta_INSS_Stage_Split_2_3 -b meta_Class -m ${confusion_path}
python "${script_path}" -i /data/gse49711 -o ${out_path} -c meta_INSS_Stage_Split_3_4 -b meta_Class -m ${confusion_path}
