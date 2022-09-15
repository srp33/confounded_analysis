#!/bin/bash

set -e

thisDir=$(dirname $0)

printf "\033[0;32mPreparing the GSE20194 dataset\033[0m\n"

#Rscript ${thisDir}/gse20194.R
Rscript ${thisDir}/create_smaller_csv.R /data/gse20194/unadjusted.csv 200 /data/gse20194/unadjusted_small.csv

#clinical_tsv="/tmp/GSE37199_Clinical.txt"
