#!/bin/bash
#SBATCH --job-name=build_base_image
#SBATCH --qos=login
#SBATCH --partition=login
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=build_base_image.out

module load apptainer

echo "Building base image with system dependencies and core R packages..."
echo "This needs to be done only once, then can be reused for fast rebuilds."

echo "Checking internet connection..."
if ! wget -q --spider http://google.com; then
    echo "ERROR: No internet connection. Apptainer build requires internet access." >&2
    exit 1
fi

# Set up directories
export APPTAINER_TMPDIR="/var/tmp/$USER/apptainer_build_tmp"
export APPTAINER_CACHEDIR="$HOME/.apptainer/cache"
mkdir -p "$APPTAINER_TMPDIR"
mkdir -p "$APPTAINER_CACHEDIR"

echo "Starting base image build at $(date)"

# Create output directory if it doesn't exist
mkdir -p ~/confounded_analysis/grp_batch_effects

# Build the base image (use --force to overwrite existing)
apptainer build --force \
    ~/confounded_analysis/grp_batch_effects/remove-batch-effects-base.sif \
    ~/confounded_analysis/apptainer/apptainer_base.def

BUILD_EXIT_CODE=$?

echo "Base image build finished at $(date) with exit code $BUILD_EXIT_CODE"

if [ $BUILD_EXIT_CODE -ne 0 ]; then
    echo "ERROR: Base image build failed with exit code $BUILD_EXIT_CODE." >&2
    exit $BUILD_EXIT_CODE
fi

echo "Base image build successful!"
echo "You can now use build_fast_image.sh for quick rebuilds."