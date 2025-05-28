#!/bin/bash

set -e

printf "\033[0;32mAdjusting the data with basic scaler\033[0m\n"

Rscript /scripts/adjust/adjust.R /data/gse20194/unadjusted.csv /data/gse20194/min_mean.csv -a min_mean -b batch
Rscript /scripts/adjust/adjust.R /data/gse24080/unadjusted.csv /data/gse24080/min_mean.csv -a min_mean -b batch
Rscript /scripts/adjust/adjust.R /data/gse49711/unadjusted.csv /data/gse49711/min_mean.csv -a min_mean -b Class