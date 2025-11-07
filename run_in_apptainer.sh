#! /bin/bash
# run_in_apptainer.sh
#
# A friendly wrapper to run commands in the Apptainer container
# under the grp_batch_effects group.
# 
# You can run it in shell mode, or run executables
# It defaults to using remove-batch-effects.sif, but other images can be specified using --image-path.
#
# It runs python, R scripts, and bash scripts automatically without specifying the executable
#
# It forwards all arguments appropriately
#
# Usage:
#   ./run_in_apptainer.sh [--image-path <path>] shell
#   ./run_in_apptainer.sh [--image-path <path>] <executable> [args...]
#   ./run_in_apptainer.sh [--image-path <path>] <script.R> [args...]
#   ./run_in_apptainer.sh [--image-path <path>] <script.py> [args...]
#   ./run_in_apptainer.sh [--image-path <path>] <script.sh> [args...]
#   ./run_in_apptainer.sh [--image-path <path>] --sbatch [sbatch-flags] <script> [script-args...]

set -euo pipefail

show_help() {
    echo "Usage:"
    echo "  $0 [--app-image-path <path>] shell"
    echo "  $0 [--app-image-path <path>] <executable> [args...]"
    echo "  $0 [--app-image-path <path>] <script.R> [args...]"
    echo "  $0 [--app-image-path <path>] <script.py> [args...]"
    echo "  $0 [--app-image-path <path>] <script.sh> [args...]"
    echo "  $0 [--app-image-path <path>] --app-sbatch [sbatch-flags] <script> [script-args...]"
    echo ""
    echo "Options:"
    echo "  --app-image-path <path>  Override the default Apptainer image path"
    echo "  --app-help               Show this help message"
    echo ""
    echo "Examples:"
    echo "  $0 --app-sbatch /scripts/evaluations/robustifying/code/3_real_data_pipe.R"
    echo "  $0 --app-sbatch --time 00:10:00 --mem 1G script.R  # Override defaults"
    echo "  $0 --app-sbatch --time 02:00:00 --cpus-per-task 4 script.py arg1 arg2"
    echo "  $0 --app-image-path /path/to/custom.sif shell"
    echo ""
    echo "Default sbatch settings: --time 1:00:00 --mem 32G --ntasks 4 --nodes 1"
}

# Parse global options manually to handle mixed option styles
TEMP_ARGS=()
CUSTOM_IMAGE_PATH=""

while [[ $# -gt 0 ]]; do
    case $1 in
        --app-image-path)
            CUSTOM_IMAGE_PATH="$2"
            shift 2
            ;;
        --app-help)
            show_help
            exit 0
            ;;
        *)
            # Save remaining arguments for later processing
            TEMP_ARGS+=("$1")
            shift
            ;;
    esac
done

# Restore remaining arguments
set -- "${TEMP_ARGS[@]}"



# Check if we have at least one command argument
if [ $# -lt 1 ]; then
    echo "Error: No command specified." >&2
    echo ""
    show_help
    exit 1
fi

MODE=$1
shift

# Load Apptainer environment (defines $APPTAINER_IMAGE)
source ~/confounded_analysis/init_apptainer.sh

# Override APPTAINER_IMAGE if custom path provided
if [[ -n "$CUSTOM_IMAGE_PATH" ]]; then
    export APPTAINER_IMAGE="$CUSTOM_IMAGE_PATH"
fi

# Helper function to detect script type and normalize path
detect_script_type() {
    if [[ $# -eq 0 || -z "${1:-}" ]]; then
        echo "Error: detect_script_type requires a script path argument" >&2
        return 1
    fi
    
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
    case "${script_path:-}" in
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
            executable="${script_path:-}"
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

    --app-sbatch)
        # Handle sbatch mode with default arguments
        declare -A SBATCH_PARAMS=(
            ["--time"]="1:00:00"
            ["--mem"]="32G"
            ["--ntasks"]="1"
            ["--cpus-per-task"]="4"
            ["--nodes"]="1"
        )
        SBATCH_FLAGS=()
        SCRIPT=""
        SCRIPT_ARGS=()
        USER_PROVIDED_JOB_NAME=false
        
        # Parse arguments to separate sbatch flags from script and its args
        while [[ $# -gt 0 ]]; do
            case $1 in
                --job-name|-J)
                    USER_PROVIDED_JOB_NAME=true
                    SBATCH_PARAMS["$1"]="$2"
                    shift 2
                    ;;
                --time|--mem|--cpus-per-task|--ntasks|--nodes|--partition|--output|-o|--error|-e|--mail-type|--mail-user|--array)
                    SBATCH_PARAMS["$1"]="$2"
                    shift 2
                    ;;
                --*)
                    # Other sbatch flags without arguments
                    SBATCH_FLAGS+=("$1")
                    shift
                    ;;
                *)
                    # First non-flag argument is the script
                    if [[ -z "$SCRIPT" ]]; then
                        SCRIPT="$1"
                    else
                        # Remaining arguments are script arguments
                        SCRIPT_ARGS+=("$1")
                    fi
                    shift
                    ;;
            esac
        done
        
        if [[ -z "$SCRIPT" ]]; then
            echo "Error: No script specified for --app-sbatch mode"
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
            SBATCH_PARAMS["--job-name"]="$SCRIPT_BASENAME"
        fi
        
        # Build final sbatch arguments array
        SBATCH_ARGS=()
        for param in "${!SBATCH_PARAMS[@]}"; do
            SBATCH_ARGS+=("$param" "${SBATCH_PARAMS[$param]}")
        done
        for flag in "${SBATCH_FLAGS[@]}"; do
            SBATCH_ARGS+=("$flag")
        done
        
        # Write sbatch script
        {
            echo '#!/bin/bash'
            echo '#SBATCH --job-name=apptainer_job'
            echo ''
            echo '# Load Apptainer environment'
            echo 'source ~/confounded_analysis/init_apptainer.sh'
            echo ''
            if [[ -n "$CUSTOM_IMAGE_PATH" ]]; then
                echo "# Override APPTAINER_IMAGE with custom path"
                echo "export APPTAINER_IMAGE=\"$CUSTOM_IMAGE_PATH\""
                echo ''
            fi
            echo '# Build command'
            echo "cmd=\"exec apptainer exec --contain \\\"\$APPTAINER_IMAGE\\\" \\\"$EXECUTABLE\\\"\""
        } > "$TEMP_SBATCH"
        
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
        if [[ -n "$EXECUTABLE" ]]; then
            detect_script_type "$EXECUTABLE"
        else
            echo "Error: No executable specified" >&2
            exit 1
        fi
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
