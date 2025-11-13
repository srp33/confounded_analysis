#!/bin/bash
# install_flexiblas.sh
# Install FlexiBLAS locally for R package compilation
# Usage: bash install_flexiblas.sh

set -e  # Exit on error

FLEXIBLAS_VERSION="3.3.1"
INSTALL_PREFIX="$HOME/.local"
BUILD_DIR="$HOME/flexiblas-build-tmp"

echo "=========================================="
echo "FlexiBLAS Local Installation"
echo "=========================================="
echo "Version: $FLEXIBLAS_VERSION"
echo "Install prefix: $INSTALL_PREFIX"
echo ""

# Check for cmake and try to load it
if ! command -v cmake &> /dev/null; then
    echo "cmake not found, attempting to load from modules..."
    
    # Try to find and load cmake module
    CMAKE_MODULE=$(module -t avail cmake 2>&1 | grep -E "^cmake/" | head -n1)
    
    if [ -n "$CMAKE_MODULE" ]; then
        echo "Loading module: $CMAKE_MODULE"
        module load "$CMAKE_MODULE"
        
        if ! command -v cmake &> /dev/null; then
            echo "ERROR: cmake still not found after loading module"
            exit 1
        fi
        echo "✓ cmake loaded successfully"
    else
        echo "ERROR: cmake not found and no cmake module available"
        echo "Check available modules: module avail cmake"
        exit 1
    fi
else
    echo "✓ cmake found: $(cmake --version | head -n1)"
fi

# Create build directory
echo "Creating build directory..."
mkdir -p "$BUILD_DIR"
cd "$BUILD_DIR"

# Download FlexiBLAS
if [ ! -f "v${FLEXIBLAS_VERSION}.tar.gz" ]; then
    echo "Downloading FlexiBLAS ${FLEXIBLAS_VERSION}..."
    wget "https://github.com/mpimd-csc/flexiblas/archive/refs/tags/v${FLEXIBLAS_VERSION}.tar.gz"
else
    echo "Using existing download: v${FLEXIBLAS_VERSION}.tar.gz"
fi

# Extract
echo "Extracting..."
tar xzf "v${FLEXIBLAS_VERSION}.tar.gz"
cd "flexiblas-${FLEXIBLAS_VERSION}"

# Configure with cmake
echo "Configuring with cmake..."
mkdir -p build
cd build
cmake -DCMAKE_INSTALL_PREFIX="$INSTALL_PREFIX" ..

# Build
echo "Building (this may take a few minutes)..."
make -j4

# Install
echo "Installing to $INSTALL_PREFIX..."
make install

# Verify installation
echo ""
echo "Verifying installation..."
if [ -f "$INSTALL_PREFIX/lib/libflexiblas.so.3" ]; then
    echo "✓ FlexiBLAS installed successfully"
    ls -lh "$INSTALL_PREFIX/lib/libflexiblas.so"*
else
    echo "ERROR: Installation verification failed"
    echo "Expected file not found: $INSTALL_PREFIX/lib/libflexiblas.so.3"
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
echo "FlexiBLAS has been installed to: $INSTALL_PREFIX"
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
echo "  3. R packages should now compile successfully"
echo ""
