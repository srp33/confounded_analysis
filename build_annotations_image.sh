#!/bin/bash
#SBATCH --job-name=build_annotations_image
#SBATCH --qos=login
#SBATCH --partition=login
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=build_annotations_image.out

module load apptainer

echo "Building annotations image from fast image..."
echo "This adds custom annotation packages from external sources."

# Check if fast image exists
FAST_IMAGE="$HOME/confounded_analysis/grp_batch_effects/remove-batch-effects-fast.sif"
if [ ! -f "$FAST_IMAGE" ]; then
    echo "ERROR: Fast image not found at $FAST_IMAGE" >&2
    echo "Please run build_fast_image.sh first to create the fast image." >&2
    exit 1
fi

echo "Checking internet connection..."
if ! wget -q --spider http://google.com; then
    echo "ERROR: No internet connection. Build requires internet access for annotation packages." >&2
    exit 1
fi

# Set up directories (minimal resources needed for annotation stage)
export APPTAINER_TMPDIR="/var/tmp/$USER/apptainer_build_tmp"
export APPTAINER_CACHEDIR="$HOME/.apptainer/cache"
mkdir -p "$APPTAINER_TMPDIR"
mkdir -p "$APPTAINER_CACHEDIR"

echo "Starting annotation image build at $(date)"

# Create output directory if it doesn't exist
mkdir -p ~/confounded_analysis/grp_batch_effects

# Build the annotation image from fast image (run from apptainer directory for relative paths)
cd ~/confounded_analysis/apptainer
apptainer build --force \
    ~/confounded_analysis/grp_batch_effects/remove-batch-effects.sif \
    apptainer_annotations.def

BUILD_EXIT_CODE=$?

echo "Annotation image build finished at $(date) with exit code $BUILD_EXIT_CODE"

# Get memory usage
echo "----------------------------------"
echo "--- PEAK MEMORY USAGE (sacct) ---"
echo "Waiting 10s for accounting database to update..."
sleep 10 
sacct -j $SLURM_JOB_ID.batch -o JobID,JobName,MaxRSS,State
echo "----------------------------------"

if [ $BUILD_EXIT_CODE -ne 0 ]; then
    echo "ERROR: Annotation image build failed with exit code $BUILD_EXIT_CODE." >&2
    echo "Check build_annotations_image.out for details." >&2
    echo "The fast image without annotations is still available at:" >&2
    echo "  $FAST_IMAGE" >&2
    exit $BUILD_EXIT_CODE
fi

echo "Annotation build successful!"
echo "Your complete container with annotations is ready at:"
echo "  ~/confounded_analysis/grp_batch_effects/remove-batch-effects.sif"
echo ""
echo "Available containers:"
echo "  Base image: ~/confounded_analysis/grp_batch_effects/remove-batch-effects-base.sif"
echo "  Expanded image: ~/confounded_analysis/grp_batch_effects/remove-batch-effects-fast.sif"
echo "  Full image with annotations: ~/confounded_analysis/grp_batch_effects/remove-batch-effects.sif"