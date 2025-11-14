# Pixi Environments

All environments are managed with [pixi](https://pixi.sh) - a fast, modern package manager for Python and R.

## Quick Start

### 1. Install pixi (one-time)
```bash
curl -fsSL https://pixi.sh/install.sh | bash
source ~/.bashrc

pixi config set --local run-post-link-scripts insecure # Might need this 
```

### 2. Set up an environment
```bash
cd environments/book_chapter

pixi install
```

### 3. Use the environment
```bash
# Enter interactive shell
pixi shell

# Or run commands directly
pixi run r
pixi run python
pixi run test  # Run tests
```

## Available Environments

### book_chapter
Full Bioconductor environment for batch correction analysis.

**Includes:**
- R 4.4 + Bioconductor packages
- Python 3.12 + scientific stack
- System libraries (OpenBLAS, GLPK, etc.)

**Usage:**
```bash
cd environments/book_chapter
pixi install
pixi shell
```

### example_env
Basic example environment with minimal dependencies.

**Includes:**
- R 4.4 + tidyverse
- Python 3.12 + numpy, pandas

**Usage:**
```bash
cd environments/example_env
pixi install
pixi shell
```

## Why Pixi?

**vs Conda:**
- ✅ 10-100x faster
- ✅ Better caching
- ✅ Simpler workflow

**vs rv + uv:**
- ✅ One tool instead of three
- ✅ Unified lock file
- ✅ Handles system libraries

**vs Conda + rv + uv:**
- ✅ Much simpler
- ✅ Faster
- ✅ Less disk space

## Common Tasks

```bash
# Install environment
pixi install

# Enter shell
pixi shell

# Run R
pixi run r

# Run Python
pixi run python

# Run tests
pixi run test

# Add a package
pixi add r-newpackage
pixi add python-package

# Update packages
pixi update

# Clean cache
pixi clean
```

## Project Structure

```
environments/
├── book_chapter/
│   ├── pixi.toml          # Dependencies
│   ├── pixi.lock          # Locked versions (auto-generated)
│   └── .pixi/             # Environment (auto-generated)
├── pixi_test/
│   └── pixi.toml
└── README.md
```

## Adding Packages

Edit `pixi.toml`:
```toml
[dependencies]
r-newpackage = "*"
python-newpackage = "*"
```

Then run:
```bash
pixi install
```

## Troubleshooting

### pixi not found
```bash
curl -fsSL https://pixi.sh/install.sh | bash
source ~/.bashrc
```

### Package not available
Check if it exists in conda-forge:
- https://anaconda.org/conda-forge

For GitHub R packages, create an install script (see `book_chapter/install_github.R`):
```r
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

pkgs <- list(
  mypackage = "user/repo"
)

for (pkg in names(pkgs)) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    remotes::install_github(pkgs[[pkg]], upgrade = "never")
  }
}
```

Then add a task to `pixi.toml`:
```toml
[tasks]
extra-install = "Rscript --vanilla install_github.R"
```

Run with:
```bash
pixi run extra-install
```

### Slow installation
First install is slow (downloads packages). Subsequent installs are fast (uses cache).

