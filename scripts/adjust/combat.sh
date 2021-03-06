#!/bin/bash

set -e

printf "\033[0;32mAdjusting the data with ComBat\033[0m\n"

Rscript /scripts/adjust/adjust.R /data/simulated_expression/unadjusted.csv /data/simulated_expression/combat.csv -b Batch
Rscript /scripts/adjust/adjust.R /data/bladderbatch/unadjusted.csv /data/bladderbatch/combat.csv -b batch
Rscript /scripts/adjust/adjust.R /data/gse37199/unadjusted.csv /data/gse37199/combat.csv -b plate
Rscript /scripts/adjust/adjust.R /data/tcga/unadjusted.csv /data/tcga/combat.csv -b CancerType
Rscript /scripts/adjust/adjust.R /data/tcga_medium/unadjusted.csv /data/tcga_medium/combat.csv -b CancerType
Rscript /scripts/adjust/adjust.R /data/tcga_small/unadjusted.csv /data/tcga_small/combat.csv -b CancerType
