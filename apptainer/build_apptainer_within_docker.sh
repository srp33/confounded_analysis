#! /bin/bash

# Exit immediately if a command exits with a non-zero status.
set -e

# Define the Docker image name.
image="build-apptainer:latest"

# Build the Docker image from the Dockerfile in the current directory.
docker build -t "$image" .

SHARED_DIR=$(pwd)

# Get the GID of the 'docker' group to ensure file permissions.
GROUP_ID=$(getent group docker | cut -d: -f3)

# Run the Docker container to build the Apptainer image.
# Use --privileged to allow container operations needed for Apptainer build
docker run --rm --privileged \
  -v "$SHARED_DIR/apptainer.def":/apptainer.def \
  -v "$SHARED_DIR/../install_packages.R":/install_packages.R \
  -v "$SHARED_DIR":/output \
  --env HOME=/ \
  --workdir / \
  "$image" \
  apptainer build --force /output/remove-batch-effects.sif /apptainer.def