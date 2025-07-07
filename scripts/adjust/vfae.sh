 #!/bin/bash

set -e

printf "\033[0;32mAdjusting the data with VFAE\033[0m\n"

python /scripts/adjust/run_vfae.py -i /data/gse20194/unadjusted.csv -o /data/gse20194/vfae.csv -b meta_batch
python /scripts/adjust/run_vfae.py -i /data/gse24080/unadjusted.csv -o /data/gse24080/vfae.csv -b meta_batch -e 4
python /scripts/adjust/run_vfae.py -i /data/gse49711/unadjusted.csv -o /data/gse49711/vfae.csv -b meta_Class -e 4
