#!/bin/bash
#SBATCH --job-name=build_combat_container
#SBATCH --qos=login
#SBATCH --partition=login
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=build_apptainer.out

module load apptainer

echo "Checking internet connection..."
if ! wget -q --spider http://google.com; then
    echo "ERROR: No internet connection. Apptainer build requires internet access." >&2
    exit 1
fi

# --- MEMORY DEBUGGING (START) ---
# Print the memory limit SLURM *thinks* it gave us (in MB)
echo "-----------------------------------"
echo "--- SLURM MEMORY LIMIT (Env Var) ---"
echo "SLURM_MEM_PER_NODE: $SLURM_MEM_PER_NODE (in MB)"
echo "-----------------------------------"

# Set Apptainer tmp directory to a disk location that is *not* on /home
# This is to avoid the 'nodev' mount option issue that caused the build to fail.
# /var/tmp is a common, safe location for this.
# We also set CACHEDIR to keep downloads persistent in /home.
export APPTAINER_TMPDIR="/var/tmp/$USER/apptainer_build_tmp"
export APPTAINER_CACHEDIR="$HOME/.apptainer/cache"
mkdir -p "$APPTAINER_TMPDIR"
mkdir -p "$APPTAINER_CACHEDIR"
echo "Set APPTAINER_TMPDIR to $APPTAINER_TMPDIR (to avoid 'nodev' on /home)"
echo "Set APPTAINER_CACHEDIR to $APPTAINER_CACHEDIR"
# --- MEMORY DEBUGGING (END) ---


echo "Starting Apptainer build at $(date)"

# Run the build in the *foreground*.
# We will use 'sacct' *after* the build finishes to get peak memory.
apptainer build \
    ~/confounded_analysis/apptainer/remove-batch-effects.sif \
    ~/confounded_analysis/apptainer/apptainer.def

BUILD_EXIT_CODE=$?

echo "Build finished at $(date) with exit code $BUILD_EXIT_CODE"

# --- PEAK MEMORY USAGE (sacct) ---
# We sleep for 10 seconds to give the SLURM accounting DB time to update
# with the memory usage from the build command that just finished.
echo "----------------------------------"
echo "--- PEAK MEMORY USAGE (sacct) ---"
echo "Fetching peak memory (MaxRSS) from SLURM accounting for Job $SLURM_JOB_ID..."
echo "Waiting 10s for accounting database to update..."
sleep 10 

echo "Final accounting (MaxRSS is peak memory usage *so far*):"
# Changed: Corrected typo. Was '.batch', now '$SLURM_JOB_ID.batch'.
sacct -j $SLURM_JOB_ID.batch -o JobID,JobName,MaxRSS,State
echo "----------------------------------"


if [ $BUILD_EXIT_CODE -ne 0 ]; then
    echo "ERROR: Apptainer build failed with exit code $BUILD_EXIT_CODE." >&2

    # Check for the wget error
    if [ $BUILD_EXIT_CODE -eq 255 ]; then
         echo "INFO: Exit code 255. This might be a wget failure." >&2
    fi

    # Check for OOM sigkill
    if [ $BUILD_EXIT_CODE -eq 137 ]; then
         echo "INFO: Exit code 137 strongly indicates an OOM SIGKILL." >&2
    fi
    exit $BUILD_EXIT_CODE
fi

echo "Build successful."

