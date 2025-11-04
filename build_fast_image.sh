#!/bin/bash
#SBATCH --job-name=build_fast_image
#SBATCH --qos=login
#SBATCH --partition=login
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=build_fast_image.out

module load apptainer

echo "Building fast image from base image..."

# Check if base image exists
BASE_IMAGE="$HOME/confounded_analysis/grp_batch_effects/remove-batch-effects-base.sif"
if [ ! -f "$BASE_IMAGE" ]; then
    echo "ERROR: Base image not found at $BASE_IMAGE" >&2
    echo "Please run build_base_image.sh first to create the base image." >&2
    exit 1
fi

echo "Checking internet connection..."
if ! wget -q --spider http://google.com; then
    echo "ERROR: No internet connection. Build requires internet access for Python packages." >&2
    exit 1
fi

# Set up directories (less memory needed for fast build)
export APPTAINER_TMPDIR="/var/tmp/$USER/apptainer_build_tmp"
export APPTAINER_CACHEDIR="$HOME/.apptainer/cache"
mkdir -p "$APPTAINER_TMPDIR"
mkdir -p "$APPTAINER_CACHEDIR"

echo "Starting fast image build at $(date)"

# Create output directory if it doesn't exist
mkdir -p ~/confounded_analysis/grp_batch_effects

# Build the final image from base (run from apptainer directory for relative paths)
cd ~/confounded_analysis/apptainer
apptainer build --force \
    ~/confounded_analysis/grp_batch_effects/remove-batch-effects-fast.sif \
    apptainer_fast.def

BUILD_EXIT_CODE=$?

echo "Fast image build finished at $(date) with exit code $BUILD_EXIT_CODE"

if [ $BUILD_EXIT_CODE -ne 0 ]; then
    echo "ERROR: Fast image build failed with exit code $BUILD_EXIT_CODE." >&2
    exit $BUILD_EXIT_CODE
fi

echo "Fast build successful!"
echo "Your container is ready at: ~/confounded_analysis/apptainer/remove-batch-effects-fast.sif"
echo "Now running build_annotations_image.sh"

sbatch build_annotations_image.sh