#! /bin/bash

set -e

image=enhanced-combination-gen

DOCKER_BUILDKIT=1 docker build -t $image .

mkdir -p data/r data/.cache data/gold data/raw_data data/combined_data data/raw_download
mkdir -p outputs/figures/classification outputs/figures/pca outputs/figures/mse_mmd 
mkdir -p outputs/metrics outputs/optimizations outputs/tables outputs/metrics/archive

#docker run -d --rm \
docker run -i -t --rm \
  --user $(id -u):$(id -g) \
  --env HOME=/ \
  -v $(pwd)/data:/data \
  -v $(pwd)/outputs:/outputs \
  -v $(pwd)/scripts:/scripts \
  -v $(pwd)/tmp:/tmp \
  -v $(pwd)/data/.cache:/.cache \
  -v $(pwd)/tmp/matplotlib/config:/.config/matplotlib \
  -v $(pwd)/tmp/matplotlib/cache:/.cache/matplotlib \
  $image \
  bash -c /scripts/all.sh

#chmod 777 outputs -R
