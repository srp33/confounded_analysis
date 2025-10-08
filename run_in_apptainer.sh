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
#   ./run_in_apptainer.sh --sbatch [sbatch-flags] <script> [script-args...]

set -euo pipefail

if [ $# -lt 1 ]; then
    echo "Usage:"
    echo "  $0 shell"
    echo "  $0 <executable> [args...]"
    echo "  $0 <script.R> [args...]"
    echo "  $0 <script.py> [args...]"
    echo "  $0 <script.sh> [args...]"
    echo "  $0 --sbatch [sbatch-flags] <script> [script-args...]"
    echo ""
    echo "Examples:"
    echo "  $0 --sbatch /scripts/evaluations/robustifying/code/3_real_data_pipe.R"
    echo "  $0 --sbatch --time 00:10:00 --mem 1G script.R  # Override defaults"
    echo "  $0 --sbatch --time 02:00:00 --cpus-per-task 4 script.py arg1 arg2"
    echo ""
    echo "Default sbatch settings: --time 1:00:00 --mem 32G --ntasks 4 --nodes 1"
    exit 1
fi

MODE=$1
shift

# Load Apptainer environment (defines $APPTAINER_IMAGE)
source ~/confounded_analysis/init_apptainer.sh

# Helper function to detect script type and normalize path
detect_script_type() {
    local script_path="$1"
    local executable=""
    local normalized_script=""
    
    # Ensure script path starts with / for apptainer container
    if [[ "$script_path" != /* ]]; then
        normalized_script="/$script_path"
    else
        normalized_script="$script_path"
    fi
    
    # Auto-detect script type
    case "$script_path" in
        *.R)
            executable="Rscript"
            ;;
        *.py)
            executable="python"
            ;;
        *.sh)
            executable="bash"
            ;;
        *)
            executable="$script_path"
            normalized_script=""
            ;;
    esac
    
    # Return values via global variables (bash limitation)
    DETECTED_EXECUTABLE="$executable"
    NORMALIZED_SCRIPT="$normalized_script"
}

case "$MODE" in
    shell|--shell)
        # Start an interactive Apptainer shell
        sg grp_batch_effects -c "apptainer shell \"$APPTAINER_IMAGE\""
        ;;

    --sbatch)
        # Handle sbatch mode with default arguments
        SBATCH_ARGS=("--time" "1:00:00" "--mem" "32G" "--ntasks" "4" "--nodes" "1")
        SCRIPT=""
        SCRIPT_ARGS=()
        USER_PROVIDED_JOB_NAME=false
        
        # Parse arguments to separate sbatch flags from script and its args
        while [[ $# -gt 0 ]]; do
            case $1 in
                --job-name|-J)
                    USER_PROVIDED_JOB_NAME=true
                    SBATCH_ARGS+=("$1")
                    if [[ $# -gt 1 && ! "$2" =~ ^-- ]]; then
                        shift
                        SBATCH_ARGS+=("$1")
                    fi
                    ;;
                --time|--mem|--cpus-per-task|--ntasks|--nodes|--partition|--output|-o|--error|-e|--mail-type|--mail-user|--array)
                    SBATCH_ARGS+=("$1")
                    if [[ $# -gt 1 && ! "$2" =~ ^-- ]]; then
                        shift
                        SBATCH_ARGS+=("$1")
                    fi
                    ;;
                --*)
                    # Other sbatch flags without arguments
                    SBATCH_ARGS+=("$1")
                    ;;
                *)
                    # First non-flag argument is the script
                    if [[ -z "$SCRIPT" ]]; then
                        SCRIPT="$1"
                    else
                        # Remaining arguments are script arguments
                        SCRIPT_ARGS+=("$1")
                    fi
                    ;;
            esac
            shift
        done
        
        if [[ -z "$SCRIPT" ]]; then
            echo "Error: No script specified for sbatch mode"
            exit 1
        fi
        
        # Create temporary sbatch script
        TEMP_SBATCH=$(mktemp /tmp/apptainer_sbatch.XXXXXX.sh)
        
        # Auto-detect script type and normalize path
        detect_script_type "$SCRIPT"
        EXECUTABLE="$DETECTED_EXECUTABLE"
        SCRIPT="$NORMALIZED_SCRIPT"
        
        # Set job name to script basename if not provided by user
        if [[ "$USER_PROVIDED_JOB_NAME" == false ]]; then
            SCRIPT_BASENAME=$(basename "$SCRIPT" | sed 's/\.[^.]*$//')  # Remove extension
            SBATCH_ARGS+=("--job-name" "$SCRIPT_BASENAME")
        fi
        
        # Write sbatch script
        cat > "$TEMP_SBATCH" << EOF
#!/bin/bash
#SBATCH --job-name=apptainer_job

# Load Apptainer environment
source ~/confounded_analysis/init_apptainer.sh

# Build command
cmd="exec apptainer exec --contain \"\$APPTAINER_IMAGE\" \"$EXECUTABLE\""
EOF
        
        if [[ -n "$SCRIPT" ]]; then
            echo "cmd=\"\$cmd \\\"$SCRIPT\\\"\"" >> "$TEMP_SBATCH"
        fi
        
        for arg in "${SCRIPT_ARGS[@]}"; do
            echo "cmd=\"\$cmd \\\"$arg\\\"\"" >> "$TEMP_SBATCH"
        done
        
        cat >> "$TEMP_SBATCH" << 'EOF'

# Execute in group context
sg grp_batch_effects -c "$cmd"
EOF
        
        # Make executable
        chmod +x "$TEMP_SBATCH"
        
        # Submit job
        echo "Submitting job with sbatch ${SBATCH_ARGS[*]} $TEMP_SBATCH"
        sbatch "${SBATCH_ARGS[@]}" "$TEMP_SBATCH"
        
        # Clean up temp file after a delay (job might not start immediately)
        (sleep 120 && rm -f "$TEMP_SBATCH") &
        ;;

    *)
        EXECUTABLE=$MODE

        # Auto-detect script type and normalize path
        detect_script_type "$EXECUTABLE"
        if [[ -n "$NORMALIZED_SCRIPT" ]]; then
            set -- "$NORMALIZED_SCRIPT" "$@"
        fi
        EXECUTABLE="$DETECTED_EXECUTABLE"

        # Build command string safely
        cmd="exec apptainer exec --contain \"$APPTAINER_IMAGE\" \"$EXECUTABLE\""
        for arg in "$@"; do
            cmd="$cmd \"$arg\""
        done

        sg grp_batch_effects -c "$cmd"
        ;;
esac
