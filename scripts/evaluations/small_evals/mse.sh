#!/bin/bash

set -e

printf "\033[0;32mCalculating MSE\033[0m\n"

out_path="/outputs/metrics/mse.csv"
pivot_path="/outputs/metrics/pivot_mse.csv"

source /scripts/evaluations/utils.sh

# Save previous file to an archive
archive_file "${out_path}"

conditional_rm "${pivot_path}"
conditional_rm "${out_path}"

script_path="$(dirname $0)/mse.py"

# Set PYTHONPATH to include the project root
export PYTHONPATH="/scripts:$PYTHONPATH"

python "${script_path}" -i /data/gold/gse49711 -o "$out_path"
python "${script_path}" -i /data/gold/gse20194 -o "$out_path"
python "${script_path}" -i /data/gold/gse24080 -o "$out_path"

python "$(dirname $0)/pivot_metics.py" -i "$out_path" -o "$pivot_path"
