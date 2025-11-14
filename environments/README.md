# Pixi Environments

All environments are managed with [pixi](https://pixi.sh) - a fast, modern package manager for Python and R.

## Quick Start

### 1. Install pixi (one-time)
```bash
curl -fsSL https://pixi.sh/install.sh | bash
source ~/.bashrc

# Read pixi documentation on this. I have a config.toml in ~/.pixi
pixi config set --local run-post-link-scripts insecure
```

### 2. Create a pixi.toml file in your project folder.
Check out the one in book_chapter or example_env.
I have one in ~/confounded_analysis/scripts/evaluations/book_chapter


### 3. Install packages. 
```bash
cd your_project_folder

pixi install
```

### 3. Install extras if needed (If packages are not in conda)
This is in my pixi.toml:
```toml
extra-install = "Rscript --vanilla install_github.R"
```

So I can run this in bash:
```
pixi run extra-install
```

### 3. Use the environment
```bash
# Enter interactive shell
pixi shell

# Or run commands directly
pixi run R
pixi run python my_python_script.py
pixi run sbatch my_sbatch_script.sh
pixi run sbatch --mem=1G --time=5m python my_python_script.py
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

pixi run python

# Run tests
pixi run <test_name>

# Run tool
pixi run <tool_name>

# Add a package—this appends it to the dependencies in the toml file
pixi add r-newpackage
pixi add python-package
pixi remove r-newpackage

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

### Package not available
Check if it exists in conda-forge or the other available channels (r, bioconda)
- https://anaconda.org/conda-forge

If not there, check again:) Then create an install script (see `book_chapter/install_github.R`):
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

