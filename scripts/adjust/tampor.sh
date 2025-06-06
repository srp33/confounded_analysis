#!/bin/bash

set -e

printf "\033[0;32mAdjusting the data with TAMPOR\033[0m\n"

Rscript /scripts/adjust/adjust.R /data/gse49711/unadjusted.csv /data/gse49711/tampor.csv -a tampor -b meta_Class &
Rscript /scripts/adjust/adjust.R /data/gse20194/unadjusted.csv /data/gse20194/tampor.csv -a tampor -b meta_batch &
Rscript /scripts/adjust/adjust.R /data/gse24080/unadjusted.csv /data/gse24080/tampor.csv -a tampor -b meta_batch &
