# Pure Conda vs Hybrid (Conda + rv + uv)

## Pure Conda Approach

**Setup:**
```bash
bash environments/conda_setup.sh book_chapter
conda activate book_chapter
```

**What you get:**
- Single conda environment with everything
- R + Python + all packages bundled
- ~2-4 GB per environment

**Pros:**
- ✅ Simple: One command to set up
- ✅ Reliable: Everything tested together
- ✅ No compilation: All pre-built binaries
- ✅ No FlexiBLAS issues: Everything bundled
- ✅ Standard conda workflow

**Cons:**
- ❌ Slow initial setup: 10-30 minutes
- ❌ Large disk usage: 2-4 GB per environment
- ❌ Slower updates: Conda dependency resolution is slow
- ❌ Limited packages: Not all R packages available in conda
- ❌ Outdated packages: Conda lags behind CRAN/Bioconductor

---

## Hybrid Approach (Current)

**Setup:**
```bash
bash environments/create_conda_r_env.sh 4.4  # One-time per R version
source environments/load_envs.sh book_chapter
```

**What you get:**
- Conda: R base + system libraries (~500 MB, shared)
- rv: R packages (~200 MB per project)
- uv: Python packages (~50 MB per project)

**Pros:**
- ✅ Fast: rv/uv are 10-100x faster than conda
- ✅ Latest packages: Direct from CRAN/PyPI/Bioconductor
- ✅ Smaller: ~750 MB per project vs 2-4 GB
- ✅ Shared R base: One conda env for multiple projects
- ✅ All packages available: Full CRAN/Bioconductor/PyPI

**Cons:**
- ❌ More complex: Three tools (conda + rv + uv)
- ❌ Some compilation: Matrix/S4Arrays compile from source
- ❌ Requires setup: Need to understand the hybrid approach

---

## Comparison Table

| Feature | Pure Conda | Hybrid (Conda + rv + uv) |
|---------|------------|--------------------------|
| **Initial setup time** | 10-30 min | 5-15 min |
| **Disk per project** | 2-4 GB | 750 MB |
| **Package updates** | Slow (conda) | Fast (rv/uv) |
| **Package availability** | Limited | Full CRAN/PyPI |
| **Package freshness** | Lags behind | Latest |
| **Compilation needed** | No | Some (Matrix, etc.) |
| **Complexity** | Simple | Moderate |
| **Shared R base** | No | Yes |
| **FlexiBLAS issues** | No | Solved (PPM+CRAN) |

---

## When to Use Each

### Use Pure Conda If:
- ✅ You want simplicity over speed
- ✅ Disk space isn't a concern
- ✅ You're comfortable with conda
- ✅ You don't need the latest packages
- ✅ You want a "batteries included" solution

### Use Hybrid If:
- ✅ You need the latest packages
- ✅ You want fast package updates
- ✅ Disk space is limited
- ✅ You have multiple projects sharing R version
- ✅ You need packages not in conda

---

## Migration Guide

### From Hybrid to Pure Conda:

1. **Create environment.yml** (already done for book_chapter)
2. **Run setup:**
   ```bash
   bash environments/conda_setup.sh book_chapter
   ```
3. **Activate:**
   ```bash
   conda activate book_chapter
   ```
4. **Remove old files** (optional):
   ```bash
   rm -rf environments/book_chapter/rv
   rm -rf environments/book_chapter/.venv
   ```

### From Pure Conda to Hybrid:

1. **Create rproject.toml and pyproject.toml**
2. **Create conda R base:**
   ```bash
   bash environments/create_conda_r_env.sh 4.4
   ```
3. **Activate:**
   ```bash
   source environments/load_envs.sh book_chapter
   ```

---

## Recommendation

**For most users: Hybrid approach**
- Faster, more flexible, smaller
- Worth the extra complexity
- Current setup works well

**For simplicity seekers: Pure conda**
- If you just want it to work
- Don't mind waiting and disk space
- Standard conda workflow

**For legacy R (< 4.3): Apptainer**
- Neither conda nor hybrid work well
- Use containerization instead
