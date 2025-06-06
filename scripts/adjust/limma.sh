#!/bin/bash

set -e

printf "\033[0;32mAdjusting the data with Limma\033[0m\n"

Rscript /scripts/adjust/adjust.R /data/gse49711/unadjusted.csv /data/gse49711/limma_target.csv -a limma -b meta_Class -c meta_INSS_Stage_Split_1_2 meta_INSS_Stage_Split_2_3 meta_INSS_Stage_Split_3_4 &
Rscript /scripts/adjust/adjust.R /data/gse20194/unadjusted.csv /data/gse20194/limma_target.csv -a limma -b meta_batch -c meta_er_status meta_her2_status meta_pr_status &
Rscript /scripts/adjust/adjust.R /data/gse24080/unadjusted.csv /data/gse24080/limma_target.csv -a limma -b meta_batch -c meta_efs_outcome_label meta_os_outcome_label &

Rscript /scripts/adjust/adjust.R /data/gse49711/unadjusted.csv /data/gse49711/limma.csv -a limma -b meta_Class &
Rscript /scripts/adjust/adjust.R /data/gse20194/unadjusted.csv /data/gse20194/limma.csv -a limma -b meta_batch &
Rscript /scripts/adjust/adjust.R /data/gse24080/unadjusted.csv /data/gse24080/limma.csv -a limma -b meta_batch &

