#! /bin/bash

set -e

image=remove-batch-effects:latest

DOCKER_BUILDKIT=1 docker build -t $image .

SHARED_DIR=$(pwd)
SCRIPTS_DIR=$(pwd)


mkdir -p $SHARED_DIR/data $SHARED_DIR/data/.cache 
mkdir -p $SHARED_DIR/outputs/figures/ $SHARED_DIR/outputs/metrics $SHARED_DIR/outputs/tables
mkdir -p $SHARED_DIR/tmp

GROUP_ID=$(getent group docker | cut -d: -f3)

#docker run -d --rm \
docker run -i -t --rm \
  --user $(id -u):$GROUP_ID \
  --env HOME=/ \
  -v $SHARED_DIR/data:/data \
  -v $SHARED_DIR/outputs:/outputs \
  -v $SCRIPTS_DIR/scripts:/scripts \
  -v $SHARED_DIR/tmp:/tmp \
  -v $SCRIPTS_DIR/data/.cache:/.cache \
  -v $SCRIPTS_DIR/tmp/matplotlib/config:/.config/matplotlib \
  -v $SCRIPTS_DIR/tmp/matplotlib/cache:/.cache/matplotlib \
  $image \
  bash -c /scripts/all.sh

#chmod 777 outputs -R