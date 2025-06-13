#!/bin/bash

set -e

printf "\033[0;32mAdjusting the data with AutoClass\033[0m\n"

# python /scripts/adjust/invert_autoclass.py

python /scripts/adjust/autoclass.py -i /data/gse49711/combat.csv -o /data/gse49711/autoclass.csv -b meta_Class
python /scripts/adjust/autoclass.py -i /data/gse20194/combat.csv -o /data/gse20194/autoclass.csv -b meta_batch
python /scripts/adjust/autoclass.py -i /data/gse24080/combat.csv -o /data/gse24080/autoclass.csv -b meta_batch

wait
