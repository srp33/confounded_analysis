#!/bin/bash
# run_with_combat_env.sh
# Wrapper for ComBat-seq environment execution
#
# This script provides backward compatibility with run_in_combat.sh
# by automatically using the combatseq.nix R environment
#
# Usage:
#   ./run_with_combat_env.sh shell
#   ./run_with_combat_env.sh <executable> [args...]
#   ./run_with_combat_env.sh <script.R> [args...]
#   ./run_with_combat_env.sh <script.py> [args...]
#   ./run_with_combat_env.sh --sbatch [sbatch-flags] <script> [script-args...]

set -euo pipefail

# Get the directory where this script is located
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Path to the main execution wrapper
RUN_WITH_ENV="${SCRIPT_DIR}/run_with_env.sh"

# Check if run_with_env.sh exists
if [ ! -f "$RUN_WITH_ENV" ]; then
    echo "Error: run_with_env.sh not found at $RUN_WITH_ENV" >&2
    exit 1
fi

# Forward all arguments to run_with_env.sh with --r-env combatseq
# This ensures the ComBat-seq environment (Bioconductor 3.11) is used
exec "$RUN_WITH_ENV" --r-env combatseq "$@"
