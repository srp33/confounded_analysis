#! /bin/bash

set -e

image=srp33/confounded-paper:version1

docker build -t $image .

mkdir -p data/simulated_expression data/bladderbatch data/gse37199 data/tcga
mkdir -p ../generated_figures ../generated_tables ../result_metrics
chmod 777 ../generated_figures ../generated_tables ../result_metrics

docker run -i -t --rm \
  --user $(id -u):$(id -g) \
  -v $(pwd)/data:/data \
  -v $(pwd)/../generated_figures:/output/figures \
  -v $(pwd)/../generated_tables:/output/tables \
  -v $(pwd)/../result_metrics:/output/metrics \
  $image
