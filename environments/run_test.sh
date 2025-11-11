#!/bin/bash
# Wrapper script to run Python environment tests

set -euo pipefail

# Source initialization
source environments/init_env.sh

# Activate Python environment
source "$PYTHON_ENV/bin/activate"

# Run the test
python environments/test_python_env.py
