#!/bin/bash

set -e

printf "\033[0;32mCalculating MMD\033[0m\n"

out_path="/outputs/metrics/mmd.csv"
pivot_path="/outputs/metrics/pivot_mmd.csv"

# Save previous files to an archive
if [ -f "${out_path}" ]; then
  mkdir -p archive
  mod_date=$(date -r "${out_path}" +%Y-%m-%d)
  filename=$(basename -- "${out_path}")
  mv "${out_path}" "archive/${mod_date}_${filename}"
fi

script_path="$(dirname $0)/mmd.py"

python "${script_path}" -i /data/gse20194 -b batch -o "$out_path"
python "${script_path}" -i /data/gse24080 -b batch -o "$out_path"
python "${script_path}" -i /data/gse49711 -b Class -o "$out_path"

python "$(dirname $0)/pivot_metics.py" -i "$out_path" -o "$pivot_path"