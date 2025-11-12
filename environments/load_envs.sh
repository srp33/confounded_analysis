#!/usr/bin/env bash
# load_envs.sh - Unified environment loader for Python (uv) and R (rv) environments
# Usage: source load_envs.sh <project_name> [options]

PROJECT_NAME="" VERBOSE=false FORCE_SYNC=false LIST_PROJECTS=false SHOW_HELP=false
HAS_UV=false HAS_RIG=false HAS_RV=false HAS_PYTHON_ENV=false HAS_R_ENV=false PYTHON_ACTIVATED=false R_ACTIVATED=false 
PROJECT_PATH="" 
ORIGINAL_DIR="$(pwd)"
trap 'cd "$ORIGINAL_DIR" 2>/dev/null || true' EXIT INT TERM

parse_arguments() {
    local positional_args=()
    while [[ $# -gt 0 ]]; do
        case $1 in
            --verbose|-v) VERBOSE=true; shift ;;
            --list|-l) LIST_PROJECTS=true; shift ;;
            --help|-h) SHOW_HELP=true; shift ;;
            -*) echo "ERROR: Unknown option: $1"; echo "Run 'source load_envs.sh --help' for usage information"; return 1 ;;
            *) positional_args+=("$1"); shift ;;
        esac
    done
    [[ ${#positional_args[@]} -gt 0 ]] && PROJECT_NAME="${positional_args[0]}"
    [[ "$VERBOSE" == true ]] && echo "[VERBOSE] PROJECT_NAME: $PROJECT_NAME, VERBOSE: $VERBOSE, FORCE_SYNC: $FORCE_SYNC"
    return 0
}

check_prerequisites() {
    [[ "$VERBOSE" == true ]] && echo "[VERBOSE] Checking for required tools..."
    command -v uv &> /dev/null && HAS_UV=true && [[ "$VERBOSE" == true ]] && echo "  ✓ uv found"
    command -v rv &> /dev/null && HAS_RV=true && [[ "$VERBOSE" == true ]] && echo "  ✓ rv found"
    # Note: R will be loaded from module system when needed, not checked here
}

detect_environments() {
    [[ "$VERBOSE" == true ]] && echo "[VERBOSE] Detecting environments for project: $PROJECT_NAME"
    PROJECT_PATH="environments/${PROJECT_NAME}"
    if [[ ! -d "$PROJECT_PATH" ]]; then
        echo -e "\nERROR: Project '$PROJECT_NAME' not found at $PROJECT_PATH/\nRun 'source load_envs.sh --list' to see available projects\n"
        return 1
    fi
    [[ -f "$PROJECT_PATH/pyproject.toml" ]] && HAS_PYTHON_ENV=true && [[ "$VERBOSE" == true ]] && echo "  ✓ Python environment detected"
    [[ -f "$PROJECT_PATH/rproject.toml" ]] && HAS_R_ENV=true && [[ "$VERBOSE" == true ]] && echo "  ✓ R environment detected"
    if [[ "$HAS_PYTHON_ENV" == false && "$HAS_R_ENV" == false ]]; then
        echo -e "\nERROR: No environments found for project '$PROJECT_NAME'\nAdd pyproject.toml or rproject.toml\n"
        return 1
    fi
}

list_available_projects() {
    echo -e "\nAvailable Projects:\n==================="
    [[ ! -d "environments" ]] && echo -e "No environments directory found.\n" && return 0
    local python_only=() r_only=() both=()
    for project_dir in environments/*/; do
        [[ ! -d "$project_dir" ]] && continue
        local project_name=$(basename "$project_dir") has_python=false has_r=false
        [[ -f "$project_dir/pyproject.toml" ]] && has_python=true
        [[ -f "$project_dir/rproject.toml" ]] && has_r=true
        [[ "$has_python" == true && "$has_r" == true ]] && both+=("$project_name") && continue
        [[ "$has_python" == true ]] && python_only+=("$project_name") && continue
        [[ "$has_r" == true ]] && r_only+=("$project_name")
    done
    [[ ${#both[@]} -gt 0 ]] && echo -e "\nPython + R:" && printf '  • %s\n' "${both[@]}"
    [[ ${#python_only[@]} -gt 0 ]] && echo -e "\nPython Only:" && printf '  • %s\n' "${python_only[@]}"
    [[ ${#r_only[@]} -gt 0 ]] && echo -e "\nR Only:" && printf '  • %s\n' "${r_only[@]}"
    local total=$((${#both[@]} + ${#python_only[@]} + ${#r_only[@]}))
    [[ $total -eq 0 ]] && echo -e "\nNo projects found." || echo -e "\nTotal: $total project(s)\nUsage: source load_envs.sh <project_name>"
    echo ""
}

deactivate_existing_environments() {
    [[ "$VERBOSE" == true ]] && echo "[VERBOSE] Checking for existing environments to deactivate..."
    
    # Deactivate Python/virtualenv
    if [[ -n "$VIRTUAL_ENV" ]]; then
        [[ "$VERBOSE" == true ]] && echo "  Deactivating existing Python environment: $VIRTUAL_ENV"
        if command -v deactivate &> /dev/null; then
            deactivate 2>/dev/null || true
        fi
        unset VIRTUAL_ENV
    fi
    
    # Deactivate Conda
    if [[ -n "$CONDA_DEFAULT_ENV" ]]; then
        [[ "$VERBOSE" == true ]] && echo "  Deactivating Conda environment: $CONDA_DEFAULT_ENV"
        if command -v conda &> /dev/null; then
            conda deactivate 2>/dev/null || true
        fi
        unset CONDA_DEFAULT_ENV
        unset CONDA_PREFIX
        unset CONDA_PYTHON_EXE
        unset CONDA_SHLVL
    fi
    
    # Clear R environment
    if [[ -n "$R_LIBS_USER" ]]; then
        [[ "$VERBOSE" == true ]] && echo "  Clearing existing R environment: $R_LIBS_USER"
        unset R_LIBS_USER
    fi
    
    [[ "$VERBOSE" == true ]] && echo "  ✓ Environment cleanup complete"
}

activate_python_env() {
    [[ "$VERBOSE" == true ]] && echo "[VERBOSE] Activating Python environment..."
    local abs_project_path="$(cd "$PROJECT_PATH" && pwd)" || { echo "ERROR: Failed to resolve $PROJECT_PATH"; return 1; }
    local needs_sync=false
    [[ ! -d "$abs_project_path/.venv" ]] && needs_sync=true && [[ "$VERBOSE" == true ]] && echo "  .venv/ missing, sync required"
    [[ "$abs_project_path/pyproject.toml" -nt "$abs_project_path/.venv" ]] && needs_sync=true && [[ "$VERBOSE" == true ]] && echo "  pyproject.toml newer"
    [[ "$needs_sync" == false ]] && echo "  ⚡ Python environment up-to-date, skipping sync"
    if [[ "$needs_sync" == true ]]; then
        echo "  Syncing Python dependencies..."
        if [[ "$VERBOSE" == true ]]; then
            (cd "$abs_project_path" && uv sync) || { echo -e "\nERROR: uv sync failed\n"; return 1; }
        else
            (cd "$abs_project_path" && uv sync --quiet 2>&1 | grep -v "^Resolved\|^Prepared\|^Installed" || true)
            [[ ${PIPESTATUS[0]} -ne 0 ]] && echo -e "\nERROR: uv sync failed\n" && return 1
        fi
    fi
    [[ ! -f "$abs_project_path/.venv/bin/activate" ]] && echo -e "\nERROR: .venv/bin/activate not found\n" && return 1
    source "$abs_project_path/.venv/bin/activate"
    [[ -z "$VIRTUAL_ENV" ]] && echo -e "\nERROR: Failed to activate Python environment\n" && return 1
    local python_version=$(python --version 2>&1)
    [[ "$VERBOSE" == true ]] && echo "  ✓ Python activated: $python_version at $VIRTUAL_ENV" || echo "  ✓ Python: $python_version"
    PYTHON_ACTIVATED=true
}

activate_r_env() {
    [[ "$VERBOSE" == true ]] && echo "[VERBOSE] Activating R environment..."
    local abs_project_path="$(cd "$PROJECT_PATH" && pwd)" || { echo "ERROR: Failed to resolve $PROJECT_PATH"; return 1; }
    
    # Extract R version from rproject.toml
    local r_version_spec=$(grep -E '^\s*r_version\s*=' "$abs_project_path/rproject.toml" | sed -E 's/.*=\s*"([^"]+)".*/\1/')
    [[ -z "$r_version_spec" ]] && r_version_spec="4.5"
    
    # Check if R is already available and at correct version
    local load_module=false
    if command -v R &> /dev/null; then
        local current_r_version=$(R --version 2>&1 | grep -oP 'R version \K[0-9]+\.[0-9]+' | head -n1)
        if [[ "$current_r_version" == "$r_version_spec"* ]]; then
            [[ "$VERBOSE" == true ]] && echo "  ⚡ R version $current_r_version already loaded, skipping module load"
        else
            [[ "$VERBOSE" == true ]] && echo "  R version $current_r_version doesn't match required $r_version_spec"
            load_module=true
        fi
    else
        [[ "$VERBOSE" == true ]] && echo "  R not found, loading from module system..."
        load_module=true
    fi
    
    # Load R from module system if needed
    if [[ "$load_module" == true ]]; then
        # Find matching R module
        local r_module=$(module -t avail r/ 2>&1 | grep -E "^r/${r_version_spec}" | head -n1)
        if [[ -z "$r_module" ]]; then
            echo -e "\nERROR: No R module matching version $r_version_spec found"
            echo "Available R modules:"
            module avail r/ 2>&1 | grep "^r/"
            return 1
        fi
        
        [[ "$VERBOSE" == true ]] && echo "  Loading module: $r_module"
        module load "$r_module" || { echo -e "\nERROR: Failed to load R module\n"; return 1; }
    fi
    
    # Check for rv directory (could be .rv or rv depending on rv version)
    local rv_dir="$abs_project_path/rv"
    [[ ! -d "$rv_dir" ]] && rv_dir="$abs_project_path/.rv"
    
    local needs_sync=false
    [[ ! -d "$rv_dir" ]] && needs_sync=true && [[ "$VERBOSE" == true ]] && echo "  rv/ missing, sync required"
    [[ "$FORCE_SYNC" == true ]] && needs_sync=true && [[ "$VERBOSE" == true ]] && echo "  --force-sync set"
    [[ -d "$rv_dir" && "$abs_project_path/rproject.toml" -nt "$rv_dir" ]] && needs_sync=true && [[ "$VERBOSE" == true ]] && echo "  rproject.toml newer"
    [[ "$needs_sync" == false ]] && echo "  ⚡ R environment up-to-date, skipping sync"
    if [[ "$needs_sync" == true ]]; then
        echo "  Syncing R dependencies..."
        if [[ "$VERBOSE" == true ]]; then
            (cd "$abs_project_path" && rv sync) || { echo -e "\nERROR: rv sync failed\n"; return 1; }
        else
            (cd "$abs_project_path" && rv sync 2>&1 | grep -E "(Installing|Error|Warning)" || true)
            [[ ${PIPESTATUS[0]} -ne 0 ]] && echo -e "\nERROR: rv sync failed\n" && return 1
        fi
    fi

    #   r: r/4.3.3-72k3zwr, r/4.4.0-ncfmhh4, r/4.4.0-o776kvt, r/4.5.0-xcvdvru, r/4.5.1-gg7txi7, r/4.5.1-hbue2wm, r/4.5.1-5sqddv2, r/4.5.1-264p7tz
    
    # Find the library path
    local rv_lib_path="$rv_dir/library"
    [[ ! -d "$rv_lib_path" ]] && rv_lib_path="$rv_dir"
    [[ ! -d "$rv_lib_path" ]] && echo -e "\nERROR: R library directory not found at $rv_lib_path\n" && return 1
    export R_LIBS_USER="$rv_lib_path"
    [[ "$VERBOSE" == true ]] && echo "  Set R_LIBS_USER=$R_LIBS_USER"
    local r_version=$(R --version 2>&1 | head -n 1)
    [[ "$VERBOSE" == true ]] && echo "  ✓ R activated: $r_version at $R_LIBS_USER" || echo "  ✓ R: $r_version"
    R_ACTIVATED=true
}

activate_environments() {
    [[ "$HAS_PYTHON_ENV" == true && "$HAS_UV" == false ]] && echo -e "\nERROR: Python env detected but uv not installed\nRun: bash /home/phr23/confounded_analysis/environments/install_managers.sh\n" && return 1
    [[ "$HAS_R_ENV" == true && "$HAS_RV" == false ]] && echo -e "\nERROR: R env detected but rv not installed\nRun: bash /home/phr23/confounded_analysis/environments/install_managers.sh\n" && return 1
    
    # Deactivate any existing environments before activating new ones
    deactivate_existing_environments
    
    echo -e "\nActivating environments for project: $PROJECT_NAME\n=================================================="
    local activation_failed=false
    [[ "$HAS_PYTHON_ENV" == true ]] && ! activate_python_env && activation_failed=true
    [[ "$HAS_R_ENV" == true ]] && ! activate_r_env && activation_failed=true
    [[ "$activation_failed" == true ]] && echo -e "\nEnvironment activation completed with errors" && return 1
    if [[ "$HAS_PYTHON_ENV" == true && "$HAS_R_ENV" == true ]]; then
        [[ "$VERBOSE" == true && -n "$VIRTUAL_ENV" && -n "$R_LIBS_USER" ]] && echo -e "\n[VERBOSE] Dual env: VIRTUAL_ENV=$VIRTUAL_ENV, R_LIBS_USER=$R_LIBS_USER"
        [[ -z "$VIRTUAL_ENV" || -z "$R_LIBS_USER" ]] && echo -e "\nWARNING: Dual environment activation incomplete"
    fi
    echo ""
}

display_status() {
    echo "Environment Status:"
    echo "==================="
    [[ "$PYTHON_ACTIVATED" == true ]] && echo "  ✓ Python environment active" && [[ -n "$VIRTUAL_ENV" ]] && echo "    Path: $VIRTUAL_ENV"
    [[ "$R_ACTIVATED" == true ]] && echo "  ✓ R environment active" && [[ -n "$R_LIBS_USER" ]] && echo "    Library: $R_LIBS_USER"
    echo -e "\nUsage:"
    if [[ "$PYTHON_ACTIVATED" == true && "$R_ACTIVATED" == true ]]; then
        echo "  • Run 'python' or 'R' - both environments active"
    elif [[ "$PYTHON_ACTIVATED" == true ]]; then
        echo "  • Run 'python' - packages from pyproject.toml available"
    elif [[ "$R_ACTIVATED" == true ]]; then
        echo "  • Run 'R' - packages from rproject.toml available"
    fi
    echo ""
}

#------------------------------------------------------------------------------
# Function: display_help
# Description: Display comprehensive help information
# Returns: 0
#------------------------------------------------------------------------------
display_help() {
    cat << 'EOF'

load_envs.sh - Unified Environment Loader
==========================================

DESCRIPTION:
    Automatically detect and activate Python (uv) and R (rv) environments
    for a given project. Supports Python-only, R-only, or dual Python+R
    configurations.

USAGE:
    source load_envs.sh <project_name> [options]

    Note: This script MUST be sourced (not executed) to modify the current
    shell environment.

ARGUMENTS:
    <project_name>      Name of the project in environments/ directory

OPTIONS:
    -h, --help          Display this help message
    -l, --list          List all available projects
    -v, --verbose       Enable verbose output showing each step

EXAMPLES:
    # Activate book_chapter environment (Python + R)
    source load_envs.sh book_chapter

    # List all available projects
    source load_envs.sh --list

    # Activate with verbose output
    source load_envs.sh book_chapter --verbose


ENVIRONMENT DETECTION:
    The script automatically detects which environments exist for a project:

    • Python environment: Detected by presence of pyproject.toml
      - Activates uv virtual environment (.venv/)
      - Sets VIRTUAL_ENV and modifies PATH

    • R environment: Detected by presence of rproject.toml
      - Activates rv library (.rv/)
      - Sets R_LIBS_USER environment variable

    • Dual environment: Both pyproject.toml and rproject.toml present
      - Activates both environments
      - Python environment activated first, then R
      - Both VIRTUAL_ENV and R_LIBS_USER are set

DIRECTORY STRUCTURE:
    environments/
    ├── book_chapter/
    │   ├── pyproject.toml      # Python dependencies (optional)
    │   ├── rproject.toml       # R dependencies (optional)
    │   ├── .venv/              # Python virtual env (created by uv)
    │   └── .rv/                # R library (created by rv)
    └── {other_projects}/

PREREQUISITES:
    • uv - Required for Python environments
      Install: curl -LsSf https://astral.sh/uv/install.sh | sh

    • rv - Required for R environments
      Install: cargo install --git https://github.com/dgkf/rv

    • R - Loaded from HPC module system
      Use: module load r/<version>

NOTES:
    • The script must be sourced, not executed:
      ✓ Correct:   source load_envs.sh book_chapter
      ✗ Incorrect: bash load_envs.sh book_chapter

    • Smart sync: Dependencies are only synced if configuration files
      are newer than the environment directories.

    • Performance: Cached environments activate in < 2 seconds.

EOF
    return 0
}

check_sourced() {
    if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
        echo -e "\nWARNING: This script must be sourced, not executed\nCorrect usage: source load_envs.sh <project_name>\nExample: source load_envs.sh book_chapter\n"
        return 1
    fi
}

main() {
    check_sourced || return 1
    parse_arguments "$@" || return 1
    [[ "$SHOW_HELP" == true ]] && display_help && return 0
    [[ "$LIST_PROJECTS" == true ]] && list_available_projects && return 0
    [[ -z "$PROJECT_NAME" ]] && echo -e "\nERROR: No project name provided\nUsage: source load_envs.sh <project_name>\nRun 'source load_envs.sh --list' to see projects\n" >&2 && return 1
    check_prerequisites
    detect_environments || return 3
    activate_environments || return 4
    display_status
}

main "$@"
