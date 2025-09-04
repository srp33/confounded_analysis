#! /bin/bash

set -e

image=remove-batch-effects:latest

DOCKER_BUILDKIT=1 docker build -t $image .

SHARED_DIR=$(pwd)/../groups/

mkdir -p $SHARED_DIR/data $SHARED_DIR/data/.cache 
mkdir -p $SHARED_DIR/outputs/figures/ $SHARED_DIR/outputs/metrics $SHARED_DIR/outputs/tables
mkdir -p $SHARED_DIR/tmp

apptainer exec \
  --bind "$(pwd)$SHARED_DIR/data":/data \
  --bind "$(pwd)$SHARED_DIR/outputs":/outputs \
  --bind "$(pwd)/scripts":/scripts \
  --bind "$(pwd)$SHARED_DIR/tmp":/tmp \
  --bind "$(pwd)/data/.cache":/.cache \
  --bind "$(pwd)/tmp/matplotlib/config":/.config/matplotlib \
  --bind "$(pwd)/tmp/matplotlib/cache":/.cache/matplotlib \
#  --env HOME=/ \
  $image \
  /scripts/all.sh
