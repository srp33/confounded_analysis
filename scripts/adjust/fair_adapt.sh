#!/bin/bash

set -e

printf "\033[0;32mAdjusting the data with Fair Adapt\033[0m\n"


python /scripts/adjust/fair_adapt.py -i /data/gse49711/unadjusted.csv -o /data/gse49711/fair_adapt.csv -b meta_Class
python /scripts/adjust/fair_adapt.py -i /data/gse20194/unadjusted.csv -o /data/gse20194/fair_adapt.csv -b meta_batch
python /scripts/adjust/fair_adapt.py -i /data/gse24080/unadjusted.csv -o /data/gse24080/fair_adapt.csv -b meta_batch
