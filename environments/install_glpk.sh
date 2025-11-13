#!/bin/bash
# install_glpk.sh
# Install GLPK (GNU Linear Programming Kit) locally for R package compilation
# Usage: bash install_glpk.sh

set -e  # Exit on error

GLPK_VERSION="5.0"
INSTALL_PREFIX="$HOME/.local"
BUILD_DIR="$HOME/glpk-build-tmp"

echo "=========================================="
echo "GLPK Local Installation"
echo "=========================================="
echo "Version: $GLPK_VERSION"
echo "Install prefix: $INSTALL_PREFIX"
echo ""

# Create build directory
echo "Creating build directory..."
mkdir -p "$BUILD_DIR"
cd "$BUILD_DIR"

# Download GLPK
if [ ! -f "glpk-${GLPK_VERSION}.tar.gz" ]; then
    echo "Downloading GLPK ${GLPK_VERSION}..."
    wget "https://ftp.gnu.org/gnu/glpk/glpk-${GLPK_VERSION}.tar.gz"
else
    echo "Using existing download: glpk-${GLPK_VERSION}.tar.gz"
fi

# Extract
echo "Extracting..."
tar xzf "glpk-${GLPK_VERSION}.tar.gz"
cd "glpk-${GLPK_VERSION}"

# Configure
echo "Configuring..."
./configure --prefix="$INSTALL_PREFIX"

# Build
echo "Building (this may take a few minutes)..."
make -j4

# Install
echo "Installing to $INSTALL_PREFIX..."
make install

# Verify installation
echo ""
echo "Verifying installation..."
if [ -f "$INSTALL_PREFIX/lib/libglpk.so" ]; then
    echo "✓ GLPK installed successfully"
    ls -lh "$INSTALL_PREFIX/lib/libglpk.so"*
    
    # Check version
    if [ -f "$INSTALL_PREFIX/bin/glpsol" ]; then
        echo ""
        "$INSTALL_PREFIX/bin/glpsol" --version
    fi
else
    echo "ERROR: Installation verification failed"
    echo "Expected file not found: $INSTALL_PREFIX/lib/libglpk.so"
    exit 1
fi

# Cleanup
echo ""
echo "Cleaning up build directory..."
cd ~
rm -rf "$BUILD_DIR"

# Display next steps
echo ""
echo "=========================================="
echo "Installation Complete!"
echo "=========================================="
echo ""
echo "GLPK has been installed to: $INSTALL_PREFIX"
echo ""
echo "The following environment variables are required:"
echo "  export LD_LIBRARY_PATH=$INSTALL_PREFIX/lib:\$LD_LIBRARY_PATH"
echo "  export PKG_CONFIG_PATH=$INSTALL_PREFIX/lib/pkgconfig:\$PKG_CONFIG_PATH"
echo ""
echo "These are already set in load_envs.sh"
echo ""
echo "Next steps:"
echo "  1. Verify paths: echo \$LD_LIBRARY_PATH"
echo "  2. Activate environment: source environments/load_envs.sh book_chapter"
echo "  3. R packages requiring GLPK should now compile successfully"
echo ""
