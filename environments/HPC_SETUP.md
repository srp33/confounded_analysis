## System Overview

- **R**: Loaded from HPC module system (I can only get R=4.4 to work)
- **uv**: User-local Python package manager
- **rv**: User-local R package manager

## Quick Start

### 1. Install package managers (one-time):
```bash
bash environments/install_managers.sh
```

### 2. Activate an environment:
```bash
source environments/load_envs.sh book_chapter
```

### 3. Deactivate:
```bash
deactivate
```

That's it! The script will:
- Load R 4.4 from modules
- Use rv to install R packages (PPM binaries + CRAN source fallback)
- Use uv to install Python packages (if pyproject.toml exists)


### Repository Configuration

```toml
repositories = [
    { alias = "PPM", url = "https://packagemanager.posit.co/cran/latest" },  # Fast binaries
    { alias = "CRAN", url = "https://cloud.r-project.org" },                 # Source fallback
    { alias = "BioCsoft", url = "https://bioconductor.org/packages/3.20/bioc" },
]
```

This gives you:
- ✅ Fast installation for most packages (PPM binaries)
- ✅ Compatibility for all packages (CRAN source fallback)
- ✅ No FlexiBLAS dependency issues

## Available Environments

### book_chapter
- **R**: 4.4 (module system)
- **Python**: 3.12 (uv)
- **Packages**: Full Bioconductor stack, tidyverse, batch correction tools

### env_template
- Like book_chapter, but only has the toml files (all you need to start)


## Usage

### Activate an environment:
```bash
source environments/load_envs.sh book_chapter
```

### List available projects:
```bash
source environments/load_envs.sh --list
```

### Verbose mode (for debugging):
```bash
source environments/load_envs.sh book_chapter --verbose
```

## How It Works

### R Module Loading

1. Script reads `r_version` from `rproject.toml` (e.g., "4.4")
2. Loads matching R module: `module load r/4.4.0-ncfmhh4`
3. R module provides:
   - R installation
   - OpenBLAS (BLAS/LAPACK)
   - System libraries in `LD_LIBRARY_PATH`

### R Package Installation (rv)

1. rv reads `rproject.toml` dependencies
2. Tries PPM first (pre-compiled binaries)
3. Falls back to CRAN (compiles from source) if PPM fails
4. Installs to `rv/library/` in project directory
5. Sets `R_LIBS_USER` to project's rv directory

## Creating New Environments

### 1. Create project directory:
```bash
mkdir environments/my_project
```

### 2. Create rproject.toml:
```toml
[project]
name = "my_project"
r_version = "4.4"

repositories = [
    { alias = "PPM", url = "https://packagemanager.posit.co/cran/latest" },
    { alias = "CRAN", url = "https://cloud.r-project.org" },
    { alias = "BioCsoft", url = "https://bioconductor.org/packages/3.20/bioc" },
]

dependencies = [
    "dplyr",
    "ggplot2",
    # ... your packages
]
```

### 3. (Optional) Create pyproject.toml for Python:
```toml
[project]
name = "my_project"
version = "0.1.0"
requires-python = ">=3.12"

dependencies = [
    "numpy",
    "pandas",
    # ... your packages
]
```

### 4. Activate:
```bash
source environments/load_envs.sh my_project
```

rv and uv will automatically install all dependencies.

## Troubleshooting

### R module not found
```
ERROR: No R module matching version X.Y found
```

**Solution:**
- Check available versions: `module avail r/`
- Update `r_version` in `rproject.toml` to match available version
- Currently supported: R 4.4+

### rv sync fails

**Check R is loaded:**
```bash
module list | grep r/
which R
```

**Try verbose mode:**
```bash
source environments/load_envs.sh book_chapter --verbose
```

**Common causes:**
- R module not loaded
- Wrong R version (need 4.4+)
- Network issues downloading packages

### Package compilation is slow

**This is normal** when using CRAN source fallback:
- PPM binaries install in seconds
- CRAN source compilation takes minutes per package
- Matrix, S4Arrays compile from source (5-10 minutes each)
- Other packages use PPM binaries (fast)

**To speed up:**
- Use multiple cores: rv uses all available cores by default
- Be patient on first sync, subsequent syncs are fast
- Consider using conda if compilation is too slow

### FlexiBLAS errors (should not happen with current setup)

If you see:
```
libflexiblas.so.3: cannot open shared object file
```

**This means:**
- You're using old configuration with PPM-only
- Update `rproject.toml` to include CRAN fallback:
```toml
repositories = [
    { alias = "PPM", url = "https://packagemanager.posit.co/cran/latest" },
    { alias = "CRAN", url = "https://cloud.r-project.org" },  # Add this
    # ... rest of repos
]
```

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
-  Why r space**: 10-50 MB per environment vs 500MB-2GB for conda
- **HPC-friendly**: Works with module system, doesn't require nested user namespaces
- **Modern tooling**: Better dependency resolution, cleaner lock files

### Trade-offs:
- **System libraries**: Must install missing libraries yourself (FlexiBLAS, GLPK, etc.)
- **Initial setup**: More upfront work to install system dependencies
- **Two tools**: Separate Python and R package managers

Once system libraries are installed, you get much faster workflows and access to the latest packages.
## Advantages:
- **Speed**: 10-100x fast