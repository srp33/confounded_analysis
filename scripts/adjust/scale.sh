#!/bin/bash

set -e

printf "\033[0;32mAdjusting the data with a linear scaler\033[0m\n"

Rscript /scripts/adjust/adjust.R /data/simulated_expression/unadjusted.csv /data/simulated_expression/scaled.csv -a scale -b Batch
Rscript /scripts/adjust/adjust.R /data/bladderbatch/unadjusted.csv /data/bladderbatch/scaled.csv -a scale -b batch
Rscript /scripts/adjust/adjust.R /data/gse37199/unadjusted.csv /data/gse37199/scaled.csv -a scale -b plate
Rscript /scripts/adjust/adjust.R /data/tcga/unadjusted.csv /data/tcga/scaled.csv -a scale -b CancerType
Rscript /scripts/adjust/adjust.R /data/tcga_medium/unadjusted.csv /data/tcga_medium/scaled.csv -a scale -b CancerType
Rscript /scripts/adjust/adjust.R /data/tcga_small/unadjusted.csv /data/tcga_small/scaled.csv -a scale -b CancerType
