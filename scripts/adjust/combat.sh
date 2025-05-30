#!/bin/bash

set -e

printf "\033[0;32mAdjusting the data with ComBat\033[0m\n"

Rscript /scripts/adjust/adjust.R /data/gse20194/unadjusted.csv /data/gse20194/combat_target.csv -a combat -b batch -c er_status her2_status pr_status
Rscript /scripts/adjust/adjust.R /data/gse49711/unadjusted.csv /data/gse49711/combat_target.csv -a combat -b Class -c efs_outcome_label os_outcome_label
Rscript /scripts/adjust/adjust.R /data/gse24080/unadjusted.csv /data/gse24080/combat_target.csv -a combat -b batch -c INSS_Stage

Rscript /scripts/adjust/adjust.R /data/gse20194/unadjusted.csv /data/gse20194/combat.csv -b batch
Rscript /scripts/adjust/adjust.R /data/gse24080/unadjusted.csv /data/gse24080/combat.csv -b batch
Rscript /scripts/adjust/adjust.R /data/gse49711/unadjusted.csv /data/gse49711/combat.csv -b Class
