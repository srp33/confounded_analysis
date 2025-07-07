#!/bin/bash

set -e

printf "\033[0;32mAdjusting the data with wasserstein based NN\033[0m\n"

python /scripts/adjust/wasserstein.py -i /data/gse20194/combat.csv -o /data/gse20194/wasserstein.csv -b meta_batch -e 100 -v
python /scripts/adjust/wasserstein.py -i /data/gse24080/combat.csv -o /data/gse24080/wasserstein.csv -b meta_batch -e 100 -v
python /scripts/adjust/wasserstein.py -i /data/gse49711/combat.csv -o /data/gse49711/wasserstein.csv -b meta_Class -e 100 -v
