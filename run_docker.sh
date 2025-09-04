#! /bin/bash

set -e

image=remove-batch-effects:latest

DOCKER_BUILDKIT=1 docker build -t $image .

SHARED_DIR=$(pwd)/../groups/

mkdir -p $SHARED_DIR/data $SHARED_DIR/data/.cache 
mkdir -p $SHARED_DIR/outputs/figures/ $SHARED_DIR/outputs/metrics $SHARED_DIR/outputs/tables
mkdir -p $SHARED_DIR/tmp

GROUP_ID=$(getent group grp_batch_effects | cut -d: -f3)

#docker run -d --rm \
docker run -i -t --rm \
  --user $(id -u):$GROUP_ID \
  --env HOME=/ \
  -v $(pwd)$SHARED_DIR/data:/data \
  -v $(pwd)$SHARED_DIR/outputs:/outputs \
  -v $(pwd)/scripts:/scripts \
  -v $(pwd)$SHARED_DIR/tmp:/tmp \
  -v $(pwd)/data/.cache:/.cache \
  -v $(pwd)/tmp/matplotlib/config:/.config/matplotlib \
  -v $(pwd)/tmp/matplotlib/cache:/.cache/matplotlib \
  $image \
  bash -c /scripts/all.sh

#chmod 777 outputs -R
