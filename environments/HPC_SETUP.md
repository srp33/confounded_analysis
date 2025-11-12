# HPC Environment Setup Guide

This guide explains how to use the environment management system on HPC clusters where you don't have sudo access.

## Key Differences from Standard Setup

- **R**: Loaded from HPC module system instead of using `rig`
- **uv**: Installed user-local for Python package management
- **rv**: Installed user-local for R package management

## Initial Setup

1. Install uv and rv (one-time):
   ```bash
   bash environments/install_managers.sh
   ```

2. Verify R is available via modules:
   ```bash
   module avail r/
   ```
   You should see versions like `r/4.5.0-xcvdvrur`, `r/4.5.1-gg7txi7`, etc.

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

### Compilation errors for specific packages
Some Bioconductor packages (like `rtracklayer`) may fail to compile due to compiler compatibility issues. This is usually not critical:
- Most packages will install successfully
- The environment is still usable
- You can try a different R version if needed
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

## Module System Commands

```bash
# List all available R modules
module avail r/

# Load a specific R version
module load r/4.5.1-gg7txi7

# Check loaded modules
module list

# Unload a module
module unload r/4.5.1-gg7txi7

# Get detailed info about a module
module spider r/4.5.1-gg7txi7
```
