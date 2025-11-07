#!/bin/bash
# run_in_combat.sh
#
# A wrapper script that runs run_in_apptainer.sh with the combat seq image
# This is equivalent to running:
#   ./run_in_apptainer.sh --app-image-path ~/groups/grp_batch_effects/combatseq_image.sif <args...>
#
# Usage:
#   ./run_in_combat.sh shell
#   ./run_in_combat.sh <executable> [args...]
#   ./run_in_combat.sh <script.R> [args...]
#   ./run_in_combat.sh <script.py> [args...]
#   ./run_in_combat.sh <script.sh> [args...]
#   ./run_in_combat.sh --app-sbatch [sbatch-flags] <script> [script-args...]

set -euo pipefail

# Get the directory where this script is located
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Execute run_in_apptainer.sh with the combat seq image path and all passed arguments
exec "$SCRIPT_DIR/run_in_apptainer.sh" --app-image-path ~/groups/grp_batch_effects/combatseq_image.sif "$@"