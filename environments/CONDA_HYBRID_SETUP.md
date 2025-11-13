# Conda + rv + uv Hybrid Setup

This setup uses the best of each tool:
- **Conda**: R versions + system libraries (FlexiBLAS, GLPK, etc.)
- **rv**: R package management (fast, up-to-date)
- **uv**: Python package management (10-100x faster than conda)

## Quick Start

### 1. (Optional) Load conda module

The script will automatically try to load conda from modules if not available. But you can load it manually:

```bash
module load miniconda3
# or
module load anaconda3
```

### 2. Create conda environments for each R version you need

For R 4.4 (book_chapter):
```bash
bash environments/create_conda_r_env.sh 4.4
```

For R 4.0 (combatseq):
```bash
bash environments/create_conda_r_env.sh 4.0
```

This creates conda environments with:
- R base installation
- System libraries (OpenBLAS, GLPK, libcurl, libxml2, etc.)
- Build tools (gfortran, make, cmake)

### 3. Activate your project environment

```bash
source environments/load_envs.sh book_chapter
```

The script will:
1. Detect R version from `rproject.toml`
2. Activate the matching conda environment (e.g., `r-4.4`)
3. Use rv to manage R packages (fast!)
4. Use uv to manage Python packages (if pyproject.toml exists)

## How It Works

### Environment Activation Flow

```
load_envs.sh
    ↓
Check rproject.toml for R version (e.g., "4.4")
    ↓
Look for conda env "rv-r4.4-syslibs"
    ↓
    ├─ Found? → conda activate rv-r4.4-syslibs
    │              ↓
    │          R + system libs from conda
    │              ↓
    │          R packages from rv (fast!)
    │
    └─ Not found? → Fall back to module system
                       ↓
                   module load r/4.4.x
```

### What Each Tool Manages

| Tool | Manages | Why |
|------|---------|-----|
| **Conda** | R base + system libraries | Bundles all dependencies, no compilation needed |
| **rv** | R packages | 10x faster than conda, latest CRAN/Bioconductor |
| **uv** | Python packages | 100x faster than conda, latest PyPI |

## Benefits

### vs Pure rv/uv
- ✅ No more missing `.so` files (OpenBLAS, GLPK, etc.)
- ✅ Easy R version switching (4.0, 4.4, 4.5)
- ✅ System libraries "just work"

### vs Pure Conda
- ✅ 10-100x faster package installation
- ✅ Latest packages from CRAN/PyPI
- ✅ Smaller disk usage (conda only for R base)
- ✅ Better dependency resolution

## Conda Environments

List created environments:
```bash
conda env list
```

You should see:
```
rv-r4.0-syslibs     /path/to/conda/envs/rv-r4.0-syslibs
rv-r4.4-syslibs     /path/to/conda/envs/rv-r4.4-syslibs
```

The naming convention `rv-r{version}-syslibs` makes it clear:
- `rv` - Used with rv package manager
- `r{version}` - R version (e.g., r4.4)
- `syslibs` - Provides R base + system libraries only

Each environment is ~500MB (much smaller than full conda R environments with all packages).

## Troubleshooting

### Conda environment not found
```
⚠ Conda environment 'rv-r4.4-syslibs' not found

Create it now? (Y/n):
```

**Solution**: 
- Press `Y` to automatically create it (recommended)
- Press `n` to fall back to module system
- Or create manually: `bash environments/create_conda_r_env.sh 4.4`

### Wrong R version
If `rproject.toml` says `r_version = "4.4.0"`, the script looks for conda env `rv-r4.4-syslibs`.

**Solution**: Either:
- Change `rproject.toml` to `r_version = "4.4"` (recommended)
- The script will auto-create the environment when you activate

### Conda not available

The script automatically tries to load conda from modules (miniconda3, anaconda3, conda).

If it still fails:
```
Conda not available, using module system
```

**Solution**: 
- Check available conda modules: `module avail conda`
- Load manually: `module load miniconda3` or `module load anaconda3`
- Or continue with module system (R from modules, no conda)

## Migration from Module-Only Setup

If you were using modules before:

1. Your existing rv environments still work
2. Create conda environments for R versions you use
3. Next time you activate, it will use conda automatically
4. You can remove locally-compiled FlexiBLAS/GLPK (conda provides them)

## Disk Space

- Conda env per R version: ~500MB
- rv packages per project: ~50-200MB
- uv packages per project: ~10-50MB

Total for 2 R versions + 2 projects: ~1.5GB

Compare to pure conda: ~4-6GB for the same setup.
