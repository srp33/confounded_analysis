#!/bin/bash
# run_with_env.sh
# Execution wrapper for uv/rix environments
#
# Replaces run_in_apptainer.sh with native environment execution
# Provides automatic script type detection and environment activation
#
# Usage:
#   ./run_with_env.sh [--r-env <path>] [--full-env] shell
#   ./run_with_env.sh [--r-env <path>] [--full-env] <executable> [args...]
#   ./run_with_env.sh [--r-env <path>] [--full-env] <script.R> [args...]
#   ./run_with_env.sh [--r-env <path>] [--full-env] <script.py> [args...]
#   ./run_with_env.sh [--r-env <path>] [--full-env] <script.sh> [args...]
#   ./run_with_env.sh [--r-env <path>] [--full-env] --sbatch [sbatch-flags] <script> [script-args...]

set -euo pipefail

# ============================================================================
# Configuration
# ============================================================================
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ANALYSIS_DIR="$(dirname "$SCRIPT_DIR")"

# Source environment initialization
if [ -f "${SCRIPT_DIR}/init_env.sh" ]; then
    source "${SCRIPT_DIR}/init_env.sh"
else
    echo "Error: init_env.sh not found at ${SCRIPT_DIR}/init_env.sh" >&2
    exit 1
fi

# Nix configuration for rix Phase 2 workflow
NIX_ROOT="/grphome/grp_batch_effects/nix"
NIX_CHROOT_CMD="$NIX_ROOT/nix-user-chroot $NIX_ROOT"

# Binary cache configuration (essential for nix-user-chroot)
NIX_CACHE_OPTS="--option substituters 'https://cache.nixos.org https://rstats-on-nix.cachix.org' --option trusted-public-keys 'cache.nixos.org-1:6NCHdD59X431o0gWypbMrAURkbJ16ZPMQFGspcDShjY= rstats-on-nix.cachix.org-1:9cJb6nqYZgKqgH5XJQN8FPkXqKlGqKqJqKqKqKqKqKo='"

# ============================================================================
# Help Message
# ============================================================================
show_help() {
    cat << EOF
Usage:
  $0 [options] shell
  $0 [options] <executable> [args...]
  $0 [options] <script.R> [args...]
  $0 [options] <script.py> [args...]
  $0 [options] <script.sh> [args...]
  $0 [options] --sbatch [sbatch-flags] <script> [script-args...]

Options:
  --r-env <name>    Use specific R environment directory (batch-effects or combatseq)
                    Default: batch-effects
  --full-env        Activate both Python and R environments
                    Use for scripts that need both (e.g., reticulate)
  --help            Show this help message

Script Type Detection:
  .py files  → Python environment (uv)
  .R files   → R environment (rix/Nix via nix-user-chroot)
  .sh files  → Both environments
  Other      → Execute directly with both environments

rix Phase 2 Workflow:
  - R execution uses nix-user-chroot wrapper for rootless Nix
  - Changes to R environment directory before running nix-shell
  - Sources Nix profile inside nix-user-chroot namespace
  - Uses explicit binary cache options for fast package downloads
  - Loads .Rprofile for library path isolation

Examples:
  # Python script (fast, Python-only)
  $0 scripts/adjust/autoclass.py

  # R script (uses rix Phase 2 with nix-user-chroot)
  $0 scripts/adjust/gmm_adjust.R

  # Mixed script (both environments, for reticulate)
  $0 --full-env scripts/mixed_analysis.R

  # Use ComBat-seq environment (directory name)
  $0 --r-env combatseq scripts/combat_analysis.R

  # SLURM job submission with R script
  $0 --sbatch --time 02:00:00 --mem 64G scripts/adjust/gmm_adjust.R

  # Interactive shell with both environments
  $0 shell

Default sbatch settings: --time 1:00:00 --mem 32G --cpus-per-task 4 --nodes 1
EOF
}

# ============================================================================
# Parse Global Options
# ============================================================================
R_ENV_NAME="batch-effects"
FULL_ENV=false
TEMP_ARGS=()

while [[ $# -gt 0 ]]; do
    case $1 in
        --r-env)
            R_ENV_NAME="$2"
            shift 2
            ;;
        --full-env)
            FULL_ENV=true
            shift
            ;;
        --help)
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

# ============================================================================
# Determine R Environment Directory (rix Phase 2 workflow)
# ============================================================================
# The --r-env flag now accepts directory names (batch-effects, combatseq)
# instead of .nix file paths. This is required for rix Phase 2 workflow
# because we need to cd into the directory before running nix-shell.
case "$R_ENV_NAME" in
    batch-effects)
        R_ENV_DIR="${ANALYSIS_DIR}/environments/r/batch-effects"
        ;;
    combatseq)
        R_ENV_DIR="${ANALYSIS_DIR}/environments/r/combatseq"
        ;;
    *)
        # Check if it's a directory path
        if [ -d "$R_ENV_NAME" ]; then
            R_ENV_DIR="$R_ENV_NAME"
        else
            echo "Error: Unknown R environment: $R_ENV_NAME" >&2
            echo "Valid options: batch-effects, combatseq, or path to R environment directory" >&2
            exit 1
        fi
        ;;
esac

# Verify R environment directory exists and has default.nix
if [ ! -d "$R_ENV_DIR" ]; then
    echo "Error: R environment directory not found: $R_ENV_DIR" >&2
    exit 1
fi

if [ ! -f "$R_ENV_DIR/default.nix" ]; then
    echo "Error: default.nix not found in $R_ENV_DIR" >&2
    echo "Run Phase 1 (authoring) first to generate default.nix" >&2
    exit 1
fi

# ============================================================================
# Environment Activation Functions
# ============================================================================

activate_python_env() {
    if [ ! -d "$PYTHON_ENV" ]; then
        echo "Error: Python environment not found at $PYTHON_ENV" >&2
        echo "Run: cd ${PYTHON_ENV_SPEC} && uv sync" >&2
        exit 1
    fi
    
    # Activate Python virtual environment
    source "$PYTHON_ENV/bin/activate"
}

activate_r_env() {
    local r_env_dir="$1"
    
    if [ ! -d "$r_env_dir" ]; then
        echo "Error: R environment directory not found: $r_env_dir" >&2
        exit 1
    fi
    
    if [ ! -f "$r_env_dir/default.nix" ]; then
        echo "Error: default.nix not found in $r_env_dir" >&2
        echo "Run Phase 1 (authoring) first to generate default.nix" >&2
        exit 1
    fi
    
    if [ ! -f "$NIX_ROOT/nix-user-chroot" ]; then
        echo "Error: nix-user-chroot not found at $NIX_ROOT/nix-user-chroot" >&2
        echo "Nix may not be installed. See environments/README.md" >&2
        exit 1
    fi
    
    # Check if pre-built result exists (nix-build optimization)
    if [ -L "$r_env_dir/result" ]; then
        # Use pre-built environment for instant activation
        NIX_SHELL_TARGET="result"
    else
        # Fall back to default.nix (slow 5-minute activation)
        NIX_SHELL_TARGET="default.nix"
        if [ "${WARN_SLOW_NIX:-1}" = "1" ]; then
            echo "Warning: Using default.nix (slow ~5min activation)" >&2
            echo "Run 'cd $r_env_dir && ../../build_nix_env.sh' to speed up" >&2
        fi
    fi
    
    # Note: Actual activation happens via nix-user-chroot + nix-shell command
    # This function just validates the environment exists and sets NIX_SHELL_TARGET
}

# ============================================================================
# Path Resolution
# ============================================================================
resolve_script_path() {
    local script_path="$1"
    
    # If path is already absolute, return it
    if [[ "$script_path" = /* ]]; then
        echo "$script_path"
        return
    fi
    
    # If path starts with ~, expand it
    if [[ "$script_path" = ~* ]]; then
        eval echo "$script_path"
        return
    fi
    
    # Otherwise, make it absolute relative to current directory
    echo "$(pwd)/$script_path"
}

# ============================================================================
# Script Type Detection
# ============================================================================
detect_script_type() {
    local script_path="$1"
    
    case "$script_path" in
        *.py)
            echo "python"
            ;;
        *.R)
            echo "r"
            ;;
        *.sh)
            echo "shell"
            ;;
        *)
            echo "unknown"
            ;;
    esac
}

# ============================================================================
# Execution Functions
# ============================================================================

execute_python_script() {
    local script="$1"
    shift
    
    activate_python_env
    python "$script" "$@"
}

execute_r_script() {
    local script="$1"
    shift
    
    activate_r_env "$R_ENV_DIR"
    
    # rix Phase 2 workflow: cd to R environment directory, then run nix-shell
    # Uses ./result if available (instant), otherwise default.nix (slow)
    $NIX_CHROOT_CMD bash -c "
        source ~/.nix-profile/etc/profile.d/nix.sh && \
        cd '$R_ENV_DIR' && \
        nix-shell '$NIX_SHELL_TARGET' $NIX_CACHE_OPTS --run 'Rscript \"$script\" $*'
    "
}

execute_shell_script() {
    local script="$1"
    shift
    
    # Shell scripts get both environments by default
    activate_python_env
    activate_r_env "$R_ENV_DIR"
    
    # rix Phase 2 workflow with Python environment available
    $NIX_CHROOT_CMD bash -c "
        source ~/.nix-profile/etc/profile.d/nix.sh && \
        cd '$R_ENV_DIR' && \
        nix-shell '$NIX_SHELL_TARGET' $NIX_CACHE_OPTS --run 'bash \"$script\" $*'
    "
}

execute_with_full_env() {
    local cmd="$1"
    shift
    
    # Activate both Python and R environments
    activate_python_env
    activate_r_env "$R_ENV_DIR"
    
    # rix Phase 2 workflow with Python environment available (for reticulate)
    $NIX_CHROOT_CMD bash -c "
        source ~/.nix-profile/etc/profile.d/nix.sh && \
        cd '$R_ENV_DIR' && \
        nix-shell '$NIX_SHELL_TARGET' $NIX_CACHE_OPTS --run 'export RETICULATE_PYTHON=\"$PYTHON_ENV/bin/python\"; $cmd $*'
    "
}

execute_direct() {
    local executable="$1"
    shift
    
    # For unknown types, activate both environments and execute directly
    activate_python_env
    activate_r_env "$R_ENV_DIR"
    
    # rix Phase 2 workflow for unknown executable types
    $NIX_CHROOT_CMD bash -c "
        source ~/.nix-profile/etc/profile.d/nix.sh && \
        cd '$R_ENV_DIR' && \
        nix-shell '$NIX_SHELL_TARGET' $NIX_CACHE_OPTS --run '$executable $*'
    "
}

# ============================================================================
# Main Execution Logic
# ============================================================================
MODE=$1
shift

case "$MODE" in
    shell|--shell)
        # Start an interactive shell with both environments
        echo "Starting interactive shell with Python and R environments..."
        echo "Python: $PYTHON_ENV"
        echo "R: $R_ENV_DIR"
        echo ""
        
        activate_python_env
        activate_r_env "$R_ENV_DIR"
        
        # rix Phase 2 workflow: interactive shell with both environments
        $NIX_CHROOT_CMD bash -c "
            source ~/.nix-profile/etc/profile.d/nix.sh && \
            cd '$R_ENV_DIR' && \
            nix-shell '$NIX_SHELL_TARGET' $NIX_CACHE_OPTS --run 'export RETICULATE_PYTHON=\"$PYTHON_ENV/bin/python\"; exec bash'
        "
        ;;

    --sbatch)
        # ====================================================================
        # SLURM Job Submission Mode
        # ====================================================================
        # Parse sbatch flags and script arguments
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
            echo "Error: No script specified for --sbatch mode" >&2
            exit 1
        fi
        
        # Resolve script path to absolute path
        SCRIPT=$(resolve_script_path "$SCRIPT")
        
        # Detect script type
        SCRIPT_TYPE=$(detect_script_type "$SCRIPT")
        
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
        
        # Create temporary sbatch script
        TEMP_SBATCH=$(mktemp /tmp/uvrix_sbatch.XXXXXX.sh)
        
        # Write sbatch script header
        cat > "$TEMP_SBATCH" << 'HEADER_EOF'
#!/bin/bash
# Auto-generated SLURM job script for uv/rix execution

set -euo pipefail

# Source environment initialization
HEADER_EOF
        
        # Add environment initialization
        echo "source ${SCRIPT_DIR}/init_env.sh" >> "$TEMP_SBATCH"
        echo "" >> "$TEMP_SBATCH"
        
        # Add R environment configuration (rix Phase 2 workflow)
        echo "# R environment configuration" >> "$TEMP_SBATCH"
        echo "R_ENV_DIR=\"$R_ENV_DIR\"" >> "$TEMP_SBATCH"
        echo "NIX_ROOT=\"$NIX_ROOT\"" >> "$TEMP_SBATCH"
        echo "NIX_CHROOT_CMD=\"\$NIX_ROOT/nix-user-chroot \$NIX_ROOT\"" >> "$TEMP_SBATCH"
        echo "NIX_CACHE_OPTS=\"$NIX_CACHE_OPTS\"" >> "$TEMP_SBATCH"
        echo "" >> "$TEMP_SBATCH"
        
        # Determine nix-shell target (result or default.nix)
        echo "# Determine nix-shell target" >> "$TEMP_SBATCH"
        echo "if [ -L \"\$R_ENV_DIR/result\" ]; then" >> "$TEMP_SBATCH"
        echo "    NIX_SHELL_TARGET=\"result\"" >> "$TEMP_SBATCH"
        echo "else" >> "$TEMP_SBATCH"
        echo "    NIX_SHELL_TARGET=\"default.nix\"" >> "$TEMP_SBATCH"
        echo "fi" >> "$TEMP_SBATCH"
        echo "" >> "$TEMP_SBATCH"
        
        # Build execution command based on script type and full-env flag
        echo "# Execute script" >> "$TEMP_SBATCH"
        
        if [ "$FULL_ENV" = true ]; then
            # Full environment: both Python and R
            cat >> "$TEMP_SBATCH" << 'FULLENV_EOF'
# Activate Python environment
source "$PYTHON_ENV/bin/activate"

# Execute in R environment with Python available
FULLENV_EOF
            
            case "$SCRIPT_TYPE" in
                python)
                    cat >> "$TEMP_SBATCH" << 'FULLENV_PYTHON_EOF'
# rix Phase 2 workflow: Python with R environment available
$NIX_CHROOT_CMD bash -c "
    source ~/.nix-profile/etc/profile.d/nix.sh && \
    cd '$R_ENV_DIR' && \
    nix-shell \"\$NIX_SHELL_TARGET\" $NIX_CACHE_OPTS --run 'export RETICULATE_PYTHON=\"$PYTHON_ENV/bin/python\"; python \"
FULLENV_PYTHON_EOF
                    echo "$SCRIPT\" ${SCRIPT_ARGS[*]}'" >> "$TEMP_SBATCH"
                    echo '"' >> "$TEMP_SBATCH"
                    ;;
                r)
                    cat >> "$TEMP_SBATCH" << 'FULLENV_R_EOF'
# rix Phase 2 workflow: R with Python environment available
$NIX_CHROOT_CMD bash -c "
    source ~/.nix-profile/etc/profile.d/nix.sh && \
    cd '$R_ENV_DIR' && \
    nix-shell \"\$NIX_SHELL_TARGET\" $NIX_CACHE_OPTS --run 'export RETICULATE_PYTHON=\"$PYTHON_ENV/bin/python\"; Rscript \"
FULLENV_R_EOF
                    echo "$SCRIPT\" ${SCRIPT_ARGS[*]}'" >> "$TEMP_SBATCH"
                    echo '"' >> "$TEMP_SBATCH"
                    ;;
                shell)
                    cat >> "$TEMP_SBATCH" << 'FULLENV_SHELL_EOF'
# rix Phase 2 workflow: Shell script with both environments
$NIX_CHROOT_CMD bash -c "
    source ~/.nix-profile/etc/profile.d/nix.sh && \
    cd '$R_ENV_DIR' && \
    nix-shell \"\$NIX_SHELL_TARGET\" $NIX_CACHE_OPTS --run 'export RETICULATE_PYTHON=\"$PYTHON_ENV/bin/python\"; bash \"
FULLENV_SHELL_EOF
                    echo "$SCRIPT\" ${SCRIPT_ARGS[*]}'" >> "$TEMP_SBATCH"
                    echo '"' >> "$TEMP_SBATCH"
                    ;;
                *)
                    cat >> "$TEMP_SBATCH" << 'FULLENV_UNKNOWN_EOF'
# rix Phase 2 workflow: Unknown type with both environments
$NIX_CHROOT_CMD bash -c "
    source ~/.nix-profile/etc/profile.d/nix.sh && \
    cd '$R_ENV_DIR' && \
    nix-shell \"\$NIX_SHELL_TARGET\" $NIX_CACHE_OPTS --run 'export RETICULATE_PYTHON=\"$PYTHON_ENV/bin/python\"; \"
FULLENV_UNKNOWN_EOF
                    echo "$SCRIPT\" ${SCRIPT_ARGS[*]}'" >> "$TEMP_SBATCH"
                    echo '"' >> "$TEMP_SBATCH"
                    ;;
            esac
        else
            # Optimized: only activate needed environment
            case "$SCRIPT_TYPE" in
                python)
                    cat >> "$TEMP_SBATCH" << 'PYTHON_EOF'
# Python-only execution (fast)
source "$PYTHON_ENV/bin/activate"
PYTHON_EOF
                    echo "python '$SCRIPT' ${SCRIPT_ARGS[*]}" >> "$TEMP_SBATCH"
                    ;;
                r)
                    cat >> "$TEMP_SBATCH" << 'R_EOF'
# R-only execution (fast, rix Phase 2 workflow)
$NIX_CHROOT_CMD bash -c "
    source ~/.nix-profile/etc/profile.d/nix.sh && \
    cd '$R_ENV_DIR' && \
    nix-shell \"\$NIX_SHELL_TARGET\" $NIX_CACHE_OPTS --run 'Rscript \"
R_EOF
                    echo "$SCRIPT\" ${SCRIPT_ARGS[*]}'" >> "$TEMP_SBATCH"
                    echo '"' >> "$TEMP_SBATCH"
                    ;;
                shell)
                    cat >> "$TEMP_SBATCH" << 'SHELL_EOF'
# Shell script execution (both environments, rix Phase 2 workflow)
source "$PYTHON_ENV/bin/activate"
$NIX_CHROOT_CMD bash -c "
    source ~/.nix-profile/etc/profile.d/nix.sh && \
    cd '$R_ENV_DIR' && \
    nix-shell \"\$NIX_SHELL_TARGET\" $NIX_CACHE_OPTS --run 'bash \"
SHELL_EOF
                    echo "$SCRIPT\" ${SCRIPT_ARGS[*]}'" >> "$TEMP_SBATCH"
                    echo '"' >> "$TEMP_SBATCH"
                    ;;
                *)
                    cat >> "$TEMP_SBATCH" << 'UNKNOWN_EOF'
# Unknown type: execute with both environments (rix Phase 2 workflow)
source "$PYTHON_ENV/bin/activate"
$NIX_CHROOT_CMD bash -c "
    source ~/.nix-profile/etc/profile.d/nix.sh && \
    cd '$R_ENV_DIR' && \
    nix-shell \"\$NIX_SHELL_TARGET\" $NIX_CACHE_OPTS --run '\"
UNKNOWN_EOF
                    echo "$SCRIPT\" ${SCRIPT_ARGS[*]}'" >> "$TEMP_SBATCH"
                    echo '"' >> "$TEMP_SBATCH"
                    ;;
            esac
        fi
        
        # Make executable
        chmod +x "$TEMP_SBATCH"
        
        # Submit job
        echo "Submitting job with sbatch ${SBATCH_ARGS[*]} $TEMP_SBATCH"
        echo "Script: $SCRIPT"
        echo "Type: $SCRIPT_TYPE"
        echo "R environment: $R_ENV_DIR"
        if [ "$FULL_ENV" = true ]; then
            echo "Mode: Full environment (Python + R, rix Phase 2)"
        else
            echo "Mode: Optimized (${SCRIPT_TYPE}-only, rix Phase 2)"
        fi
        echo ""
        
        sbatch "${SBATCH_ARGS[@]}" "$TEMP_SBATCH"
        
        # Clean up temp file after a delay (job might not start immediately)
        (sleep 120 && rm -f "$TEMP_SBATCH") &
        ;;

    *)
        # Resolve script path to absolute path
        SCRIPT_PATH=$(resolve_script_path "$MODE")
        
        # Detect script type and execute appropriately
        SCRIPT_TYPE=$(detect_script_type "$SCRIPT_PATH")
        
        if [ "$FULL_ENV" = true ]; then
            # Force full environment activation
            case "$SCRIPT_TYPE" in
                python)
                    execute_with_full_env "python '$SCRIPT_PATH'" "$@"
                    ;;
                r)
                    execute_with_full_env "Rscript '$SCRIPT_PATH'" "$@"
                    ;;
                shell)
                    execute_with_full_env "bash '$SCRIPT_PATH'" "$@"
                    ;;
                *)
                    execute_with_full_env "'$SCRIPT_PATH'" "$@"
                    ;;
            esac
        else
            # Optimize: only activate needed environment
            case "$SCRIPT_TYPE" in
                python)
                    execute_python_script "$SCRIPT_PATH" "$@"
                    ;;
                r)
                    execute_r_script "$SCRIPT_PATH" "$@"
                    ;;
                shell)
                    execute_shell_script "$SCRIPT_PATH" "$@"
                    ;;
                *)
                    # Unknown type: try direct execution with both envs
                    execute_direct "$SCRIPT_PATH" "$@"
                    ;;
            esac
        fi
        ;;
esac
