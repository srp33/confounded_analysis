#!/bin/bash

set -e

printf "\033[0;32mDownloading and tidying the bladderbatch dataset\033[0m\n"

Rscript bladderbatch.R /data/bladderbatch/unadjusted.csv
