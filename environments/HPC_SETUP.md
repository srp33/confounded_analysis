## System Overview

- **R**: Loaded from HPC module system (I can only get R=4.4 to work)
- **uv**: User-local Python package manager
- **rv**: User-local R package manager

## Quick Start

### 1. Install package managers (one-time):
```bash
bash environments/install_managers.sh
```

### 2. Edit Toml files
See the env_template folder for a base, see book_chapter for example

### 3. Activate an environment:
```bash
source environments/load_envs.sh book_chapter
```

### 4. Run file using environment:


### 5. Deactivate:
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
- **Initial setup**: More upfront work to install system dependencies
- **Two tools**: Separate Python and R package managers

Once system libraries are installed, you get much faster workflows and access to the latest packages.
