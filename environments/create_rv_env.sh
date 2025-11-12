#!/bin/bash
# create_rv_env.sh
# General script to create and update a self-contained R project using module R and rv.
# Usage: bash create_rv_env.sh [project_directory]

# --- Parse arguments ---
PROJECT_DIR="${1:-.}"

# --- Navigate to project directory ---
if [ ! -d "$PROJECT_DIR" ]; then
  echo "ERROR: Directory '$PROJECT_DIR' does not exist"
  exit 1
fi

cd "$PROJECT_DIR" || exit 1
PROJECT_DIR="$(pwd)"  # Get absolute path

echo "Setting up R environment in: $PROJECT_DIR"

# --- Ensure we are in a valid rv project ---
if [ ! -f "rproject.toml" ]; then
  echo "ERROR: rproject.toml not found in $PROJECT_DIR"
  echo "This script requires a project directory containing rproject.toml"
  exit 1
fi

# --- Extract R version from rproject.toml ---
R_VERSION=$(grep -E '^\s*r_version\s*=' "rproject.toml" | sed -E 's/.*=\s*"([^"]+)".*/\1/')

if [ -z "$R_VERSION" ]; then
  echo "WARNING: r_version not found in rproject.toml, defaulting to 4.5"
  R_VERSION="4.5"
else
  echo "Required R version from rproject.toml: $R_VERSION"
fi

# --- Check if R is already available and at correct version ---
LOAD_MODULE=false

if command -v R &> /dev/null; then
  # R is available, check version
  CURRENT_R_VERSION=$(R --version 2>&1 | grep -oP 'R version \K[0-9]+\.[0-9]+' | head -n1)
  echo "Found R version $CURRENT_R_VERSION"
  
  if [[ "$CURRENT_R_VERSION" == "$R_VERSION"* ]]; then
    echo "R version $CURRENT_R_VERSION matches required version $R_VERSION"
  else
    echo "R version $CURRENT_R_VERSION does not match required version $R_VERSION"
    LOAD_MODULE=true
  fi
else
  echo "R not found in PATH"
  LOAD_MODULE=true
fi

# --- Load R from module system if needed ---
if [ "$LOAD_MODULE" = true ]; then
  echo "--- Loading R version $R_VERSION from module system ---"
  # Find the best matching R module
  R_MODULE=$(module -t avail r/ 2>&1 | grep -E "^r/${R_VERSION}" | head -n1)

  if [ -z "$R_MODULE" ]; then
    echo "ERROR: No R module matching version $R_VERSION found"
    echo "Available R modules:"
    module avail r/ 2>&1 | grep "^r/"
    exit 1
  fi

  echo "Loading module: $R_MODULE"
  module load "$R_MODULE"

  # Verify R is available
  if ! command -v R &> /dev/null; then
    echo "ERROR: R not found after loading module"
    exit 1
  fi

  echo "R loaded successfully: $(R --version | head -n1)"
fi

# --- Check if rv is available ---
if ! command -v rv &> /dev/null; then
  echo "ERROR: rv not found in PATH"
  echo "Install rv first: cargo install --git https://github.com/dgkf/rv"
  echo "Or run: bash environments/install_managers.sh"
  exit 1
fi

# --- Initialize or update rv environment ---
if [ ! -d ".rv" ]; then
  echo "--- No rv environment detected. Initializing rv project ---"
  rv init .
else
  echo "--- Existing rv environment detected ---"
fi

# --- Sync the Environment ---
echo -e "\n--- Planning environment (dry run) ---"
rv plan || echo "rv plan failed (check rproject.toml configuration)"

echo -e "\n--- Syncing environment (installing packages) ---"
echo "This may take several minutes for the first sync..."
rv sync || echo "rv sync failed (check that rproject.toml is valid)"

echo -e "\n✅ Done. R environment ready at $PROJECT_DIR/.rv"
