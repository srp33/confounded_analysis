#! /bin/bash
# run_in_apptainer.sh
#
# A friendly wrapper to run commands in the Apptainer container
# under the grp_batch_effects group.
#
# Usage:
#   ./run_in_apptainer.sh shell
#   ./run_in_apptainer.sh <executable> [args...]
#   ./run_in_apptainer.sh <script.R> [args...]
#   ./run_in_apptainer.sh <script.py> [args...]
#   ./run_in_apptainer.sh <script.sh> [args...]

set -euo pipefail

if [ $# -lt 1 ]; then
    echo "Usage:"
    echo "  $0 shell"
    echo "  $0 <executable> [args...]"
    echo "  $0 <script.R> [args...]"
    echo "  $0 <script.py> [args...]"
    echo "  $0 <script.sh> [args...]"
    exit 1
fi

MODE=$1
shift

# Load Apptainer environment (defines $APPTAINER_IMAGE)
source ~/confounded_analysis/init_apptainer.sh

case "$MODE" in
    shell)
        # Start an interactive Apptainer shell
        declare -a CMD=( "apptainer" "shell" "$APPTAINER_IMAGE" )
        sg grp_batch_effects "${CMD[@]}"
        ;;

    *)
        EXECUTABLE=$MODE

        # Auto-detect common script extensions
        if [[ "$EXECUTABLE" == *.R ]]; then
            set -- "$EXECUTABLE" "$@"
            EXECUTABLE="Rscript"
        elif [[ "$EXECUTABLE" == *.py ]]; then
            set -- "$EXECUTABLE" "$@"
            EXECUTABLE="python"
        elif [[ "$EXECUTABLE" == *.sh ]]; then
            set -- "$EXECUTABLE" "$@"
            EXECUTABLE="bash"
        fi

        # Build command string safely
        cmd="exec apptainer exec --contain \"$APPTAINER_IMAGE\" \"$EXECUTABLE\""
        for arg in "$@"; do
            cmd="$cmd \"$arg\""
        done

        sg grp_batch_effects -c "$cmd"
        ;;
esac
