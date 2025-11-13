# HPC Environment Setup Guide

This guide explains how to use the environment management system on HPC clusters where you don't have sudo access.

## Key Differences from Standard Setup

- **R**: Loaded from HPC module system instead of using `rig`
- **uv**: Installed user-local for Python package management
- **rv**: Installed user-local for R package management
- **System libraries**: Installed to `~/.local` when not available via modules

## Initial Setup

### 1. Install package managers (one-time):
```bash
bash environments/install_managers.sh
```

### 2. Install required system libraries (one-time):

Some R packages require system libraries that may not be available on your HPC system. Install them locally:

```bash
# FlexiBLAS - Required for Matrix, S4Arrays, and other linear algebra packages
bash environments/install_flexiblas.sh

# GLPK - Required for igraph (used by batchelor and network analysis packages)
bash environments/install_glpk.sh
```

These libraries will be installed to `~/.local/lib` and automatically found by `load_envs.sh`.

### 3. Verify R is available via modules:
```bash
module avail r/
```
You should see versions like `r/4.4.0-ncfmhh4`, `r/4.5.1-gg7txi7`, etc.

## Using Environments

### Activate an environment:
```bash
source environments/load_envs.sh book_chapter
```

The script will:
1. Automatically load the appropriate R module based on `rproject.toml`
2. Activate the Python virtual environment (if present)
3. Set up the R library path for project-specific packages

### List available projects:
```bash
source environments/load_envs.sh --list
```

### Force dependency sync:
```bash
source environments/load_envs.sh book_chapter --force-sync
```

## How R Module Loading Works

When you activate an R environment:

1. The script reads the `r_version` from `rproject.toml` (e.g., "4.5")
2. It searches for matching modules: `module avail r/4.5*`
3. It loads the first matching module
4. R packages are installed to `.rv/` in your project directory

## Creating New R Environments

```bash
# From anywhere, specify the project directory
bash environments/create_rv_env.sh environments/book_chapter

# Or from environments
bash create_rv_env.sh book_chapter
```

This will:
1. Read the R version from `rproject.toml`
2. Load R from the module system (if not already loaded)
3. Initialize an rv project
4. Install packages listed in `rproject.toml`

## Troubleshooting

### R module not found
If you get "No R module matching version X.Y found":
- Check available versions: `module avail r/`
- Update `r_version` in `rproject.toml` to match an available version

### rv sync fails
- Ensure R is loaded: `module list | grep r/`
- Check R is accessible: `which R`
- Try verbose mode: `source environments/load_envs.sh book_chapter --verbose`

### Missing system library errors

If you see errors like:
```
unable to load shared object: libflexiblas.so.3: cannot open shared object file
unable to load shared object: libglpk.so.40: cannot open shared object file
```

**Solution**: Install the missing library locally:
- For FlexiBLAS: `bash environments/install_flexiblas.sh`
- For GLPK: `bash environments/install_glpk.sh`

**Why this happens**: R packages compiled from source link against system libraries. If those libraries aren't available via the module system, you need to install them locally.

**Common missing libraries**:
- `libflexiblas.so.3` - Required by: Matrix, S4Arrays, many linear algebra packages
- `libglpk.so.40` - Required by: igraph (used by batchelor, network analysis)

### Compilation errors for specific packages
Some Bioconductor packages may fail to compile due to missing dependencies:
- Check the error message for missing `.so` files
- Install the required system library (see above)
- If the library isn't critical, comment out the package in `rproject.toml`
- Failed packages can be installed manually later if required

## Performance Optimization for Many Jobs

The module system can be slow (10+ seconds per load). Here are strategies to avoid repeated module loading:

### 1. Pre-load R in your job script or interactive session
```bash
# Load R once at the start of your session
module load r/4.5.1-gg7txi7

# Then activate environments - they'll skip module loading
source environments/load_envs.sh book_chapter
```

### 2. For Snakemake/workflow systems
Load R once in the main script before launching the workflow:
```bash
#!/bin/bash
module load r/4.5.1-gg7txi7  # Load once
source environments/load_envs.sh book_chapter
snakemake --cores 8  # All jobs inherit the loaded module
```

### 3. For SLURM job arrays
Include module load in the SLURM script header:
```bash
#!/bin/bash
#SBATCH --array=1-100
#SBATCH --time=1:00:00

module load r/4.5.1-gg7txi7  # Loaded once per job
source environments/load_envs.sh book_chapter

# Your analysis here
```

### 4. Cache the module path
If you frequently use the same R version, you can set it permanently:
```bash
# Add to ~/.bashrc or ~/.bash_profile
module load r/4.5.1-gg7txi7
```

### How the optimization works
The updated `load_envs.sh` and `create_rv_env.sh` scripts now:
- Check if R is already available before attempting module load
- Verify the R version matches requirements
- Skip module loading entirely if R is already correct (saves 10+ seconds)

## Why rv + uv instead of Conda?

### Advantages:
- **Speed**: 10-100x faster package installation and environment activation
- **Native ecosystems**: Direct access to PyPI (500k+ packages) and CRAN/Bioconductor (20k+ packages)
- **Disk space**: 10-50 MB per environment vs 500MB-2GB for conda
- **HPC-friendly**: Works with module system, doesn't require nested user namespaces
- **Modern tooling**: Better dependency resolution, cleaner lock files

### Trade-offs:
- **System libraries**: Must install missing libraries yourself (FlexiBLAS, GLPK, etc.)
- **Initial setup**: More upfront work to install system dependencies
- **Two tools**: Separate Python and R package managers

Once system libraries are installed, you get much faster workflows and access to the latest packages.

## Module System Commands

```bash
# List all available R modules
module avail r/

# Load a specific R version
module load r/4.4.0-ncfmhh4

# Check loaded modules
module list

# Unload a module
module unload r/4.4.0-ncfmhh4

# Get detailed info about a module
module spider r/4.4.0-ncfmhh4
```

## Installed System Libraries

The following libraries are installed to `~/.local` and automatically found via `LD_LIBRARY_PATH`:

- **FlexiBLAS** (`libflexiblas.so.3`) - BLAS/LAPACK wrapper for linear algebra
- **GLPK** (`libglpk.so.40`) - GNU Linear Programming Kit for optimization

To verify installation:
```bash
ls ~/.local/lib/libflexiblas* ~/.local/lib/libglpk*
echo $LD_LIBRARY_PATH  # Should include ~/.local/lib
```
