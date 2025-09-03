#!/bin/bash

set -e

printf "\033[0;32mCalculating MMD\033[0m\n"

out_path="/outputs/metrics/mmd.csv"
pivot_path="/outputs/metrics/pivot_mmd.csv"

source /scripts/metrics/utils.sh

# Save previous file to an archive
archive_file "${out_path}"

conditional_rm "${pivot_path}"
conditional_rm "${out_path}"

script_path="$(dirname $0)/mmd.py"

python "${script_path}" -i /data/gold/gse49711 -b meta_Sex -o "$out_path"
python "${script_path}" -i /data/gold/gse20194 -b meta_batch -o "$out_path"
python "${script_path}" -i /data/gold/gse24080 -b meta_batch -o "$out_path"

python "$(dirname $0)/pivot_metics.py" -i "$out_path" -o "$pivot_path"