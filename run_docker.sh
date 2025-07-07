#! /bin/bash

set -e

image=srp33/confounded-paper:version1

docker build -t $image .

mkdir -p data/gse20194 data/gse24080 data/gse49711 data/r data/simple2d
mkdir -p outputs/figures outputs/metrics outputs/optimizations outputs/tables outputs/metrics/archive

#docker run -d --rm \
docker run -i -t --rm \
  --user $(id -u):$(id -g) \
  --env HOME=/ \
  -v $(pwd)/data:/data \
  -v $(pwd)/outputs:/outputs \
  -v $(pwd)/scripts:/scripts \
  -v $(pwd)/tmp:/tmp \
  -v $(pwd)/tmp/matplotlib/config:/.config/matplotlib \
  -v $(pwd)/tmp/matplotlib/cache:/.cache/matplotlib \
  $image \
  bash -c /scripts/all.sh

#chmod 777 outputs -R
