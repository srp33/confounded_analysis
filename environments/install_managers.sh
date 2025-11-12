#!/bin/bash
# Unified environment bootstrap script for user-local installation of uv and rv.
# For HPC environments where R is available via module system.

set -euo pipefail

echo "--- Unified Environment Installer ---"
echo "(uv, rv; user-local, no sudo)"
echo "Note: R will be loaded from module system"

# --- PATH setup ---
case ":$PATH:" in
  *":$HOME/.local/bin:"*) ;;
  *) export PATH="$HOME/.local/bin:$PATH" ;;
esac

mkdir -p "$HOME/.local/bin" "$HOME/.local/share" "$HOME/.local/src"

# --- Check uv installation ---
if ! command -v uv &> /dev/null; then
    echo "Installing uv..."
    curl -LsSf https://astral.sh/uv/install.sh | sh
else
    echo "✅ uv is already installed: $(uv --version)"
fi


# --- Install rv (R package manager) ---
if ! command -v rv &> /dev/null; then
    echo "Installing rv..."
    curl -sSL https://raw.githubusercontent.com/A2-ai/rv/refs/heads/main/scripts/install.sh | bash
else
    echo "✅ rv is already installed: $(rv --version)"
fi


# --- Verify installations ---
echo
echo "--- Verifying installations ---"
uv --version || echo "uv missing"
rv --version || echo "rv missing"

echo
echo "--- Environment setup complete ---"
echo "Note: Use 'module load r/<version>' to load R from your HPC module system"
