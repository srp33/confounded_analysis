#!/bin/bash

set -e

printf "\033[0;32mCalculating Single Metric Comparisons\033[0m\n"


python single_metric.py -i "/output/metrics/" -o "/output/metrics/"