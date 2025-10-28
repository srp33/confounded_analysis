#!/bin/bash
#SBATCH --job-name=build_base_image
#SBATCH --qos=login
#SBATCH --partition=login
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=build_base_image.out

module load apptainer

LOGDIR="./build_logs"
DIAG_LOG="$LOGDIR/build_diagnostics.log"
mkdir -p "$LOGDIR"

echo "=== Starting diagnostic logging ==="
echo "Timestamp: $(date)" > "$DIAG_LOG"
echo "Hostname: $(hostname)" >> "$DIAG_LOG"
echo "Node resources:" >> "$DIAG_LOG"
free -h >> "$DIAG_LOG"
df -h >> "$DIAG_LOG"
echo "CPUs: $(nproc)" >> "$DIAG_LOG"
echo "" >> "$DIAG_LOG"

echo "Building base image with system dependencies and core R packages..."
echo "Checking DNS resolution..."
if ! getent hosts google.com > /dev/null; then
    echo "DNS ERROR: Cannot resolve google.com" | tee -a "$DIAG_LOG" >&2
fi

echo "Checking internet connectivity..."
if ! wget -q --spider https://cloud.r-project.org; then
    echo "NETWORK ERROR: Cannot reach CRAN or internet may be blocked." | tee -a "$DIAG_LOG" >&2
fi

# Set up directories
export APPTAINER_TMPDIR="/var/tmp/$USER/apptainer_build_tmp"
export APPTAINER_CACHEDIR="$HOME/.apptainer/cache"
mkdir -p "$APPTAINER_TMPDIR" "$APPTAINER_CACHEDIR"

echo "Starting base image build at $(date)"

TARGET_SIF=~/confounded_analysis/grp_batch_effects/remove-batch-effects-base.sif
DEF_FILE=~/confounded_analysis/apptainer/apptainer_base.def

mkdir -p "$(dirname "$TARGET_SIF")"

apptainer build --force "$TARGET_SIF" "$DEF_FILE"
BUILD_EXIT_CODE=$?

echo "Base image build finished at $(date) with exit code $BUILD_EXIT_CODE"

if [ $BUILD_EXIT_CODE -ne 0 ]; then
    echo "ERROR: Base image build failed." | tee -a "$DIAG_LOG" >&2

    echo "" >> "$DIAG_LOG"
    echo "=== Collecting diagnostic information ===" | tee -a "$DIAG_LOG"
    echo "" >> "$DIAG_LOG"

    echo "--- dmesg OOM scan ---" >> "$DIAG_LOG"
    dmesg | tail -200 | grep -i -E "killed process|out of memory|oom" >> "$DIAG_LOG" || echo "No OOM messages detected" >> "$DIAG_LOG"

    echo "" >> "$DIAG_LOG"
    echo "--- Checking if Apptainer was killed by SIGKILL (likely OOM) ---" >> "$DIAG_LOG"
    if [ $BUILD_EXIT_CODE -eq 137 ]; then
        echo "Apptainer terminated by SIGKILL (OOM likely)" >> "$DIAG_LOG"
    fi

    echo "" >> "$DIAG_LOG"
    echo "--- Recent syslog snapshot (if accessible) ---" >> "$DIAG_LOG"
    if [ -r /var/log/syslog ]; then
        tail -200 /var/log/syslog >> "$DIAG_LOG"
    else
        echo "Syslog not readable" >> "$DIAG_LOG"
    fi

    echo "" >> "$DIAG_LOG"
    echo "--- Network route information ---" >> "$DIAG_LOG"
    ip route >> "$DIAG_LOG"

    echo "" >> "$DIAG_LOG"
    echo "--- DNS Nameservers ---" >> "$DIAG_LOG"
    cat /etc/resolv.conf >> "$DIAG_LOG"

    echo "" >> "$DIAG_LOG"
    echo "Creating diagnostic archive..."
    tar -czf build_failure_diagnostics.tgz "$LOGDIR"

    echo "All diagnostics captured in build_failure_diagnostics.tgz" >&2
    exit $BUILD_EXIT_CODE
fi

echo "Base image build successful!"
echo "Logs written to $DIAG_LOG"
echo "You can now use build_fast_image.sh for quick rebuilds."
