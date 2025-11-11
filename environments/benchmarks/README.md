# Performance Benchmarks

This directory contains scripts for benchmarking the performance of Apptainer vs uv/rix environments for the batch effect correction analysis pipeline.

## Overview

The benchmarks measure:

1. **Job Startup Times** - Time to activate environments and start executing code
2. **Package Import Times** - Time to import/load commonly used packages
3. **SLURM Array Jobs** - Throughput and reliability of array job execution

## Quick Start

Run all benchmarks:

```bash
cd environments/benchmarks
bash run_all_benchmarks.sh
```

This will:
- Run startup time benchmarks (10 iterations each)
- Run import time benchmarks on login node (10 iterations each)
- Submit import time benchmarks to compute node via SLURM
- Submit array job benchmarks (10 and 100 task arrays)
- Generate a comprehensive summary report

Results will be saved to `results_YYYYMMDD_HHMMSS/` with a summary in `BENCHMARK_SUMMARY.md`.

## Individual Benchmarks

### 1. Job Startup Times

Measures environment activation overhead:

```bash
bash benchmark_startup.sh
```

**Measures:**
- Apptainer container startup (with bind mounts)
- uv Python environment activation
- nix-shell R environment activation
- Combined Python + R activation

**Output:** `startup_times.csv`

### 2. Package Import Times

Measures package import/load performance:

```bash
# On login node
bash benchmark_imports.sh login

# On compute node (via SLURM)
sbatch --wrap="bash benchmark_imports.sh compute"
```

**Measures:**
- Python: NumPy, PyTorch
- R: tidyverse, limma (Bioconductor)

**Output:** `import_times.csv`

### 3. SLURM Array Jobs

Tests array job execution:

```bash
bash benchmark_array_jobs.sh
```

**Tests:**
- Small array (10 tasks) with Python and R
- Medium array (100 tasks) with Python and R
- Compares Apptainer vs uv/rix

**Output:** `array_job_results/` directory with:
- Job statistics for each test
- Individual task outputs
- Summary report

## Requirements

### For All Benchmarks

- Access to BYU RC cluster with SLURM
- Apptainer image at `~/groups/grp_batch_effects/remove-batch-effects.sif`
- Python environment set up (Task 7 complete)
- `bc` command for statistics calculations

### For R Benchmarks

- R environment built (Task 8.3 complete)
- Nix installed at `/grphome/grp_batch_effects/nix/`
- R environment directory at `environments/r/batch-effects/`

If R environment is not available, R benchmarks will be skipped automatically.

## Interpreting Results

### Startup Times

- **Lower is better** - Faster environment activation means less overhead
- **Target:** < 500ms for production use
- **Comparison:** Apptainer typically 1-2s, uv/rix typically 100-500ms

### Import Times

- **Lower is better** - Faster imports mean quicker job execution
- **Comparison:** Should be comparable between Apptainer and uv/rix
- **Note:** First import may be slower due to caching

### Array Jobs

- **Success rate** - All tasks should complete successfully
- **Throughput** - Total wall time from first task start to last task end
- **Comparison:** Both systems should handle array jobs reliably

## Expected Results

Based on design estimates:

| Metric | Apptainer | uv/rix | Improvement |
|--------|-----------|--------|-------------|
| Job startup | 1-2s | 100-500ms | 2-10x faster |
| Python imports | Baseline | Comparable | Similar |
| R package loads | Baseline | Comparable | Similar |
| Array job throughput | Baseline | Comparable | Similar |
| Storage | 41 GB | 18-27 GB | 40-50% reduction |

## Troubleshooting

### "R environment not found"

R benchmarks require the R environment to be built first:

```bash
cd environments/r/batch-effects
bash run_generator.sh  # Phase 1: Generate default.nix
bash build_with_cache.sh  # Phase 2: Build environment
```

### "Nix not found"

Ensure Nix is installed:

```bash
ls -la /grphome/grp_batch_effects/nix/nix-user-chroot
```

If not found, complete Task 1 (Nix installation).

### "Python environment not found"

Ensure Python environment is set up:

```bash
cd environments/python
uv sync
```

### SLURM jobs fail

Check SLURM queue and job logs:

```bash
squeue -u $USER
sacct -j <job_id>
cat <output_file>
```

## Files

- `benchmark_startup.sh` - Job startup time benchmarks
- `benchmark_imports.sh` - Package import time benchmarks
- `benchmark_array_jobs.sh` - SLURM array job benchmarks
- `run_all_benchmarks.sh` - Master script to run all benchmarks
- `README.md` - This file

## Results

Results are saved to timestamped directories:

```
results_YYYYMMDD_HHMMSS/
├── BENCHMARK_SUMMARY.md          # Summary report
├── startup_times.csv             # Startup time data
├── import_times_login.csv        # Import times (login node)
├── import_times_combined.csv     # Import times (all nodes)
├── array_job_results/            # Array job results
│   ├── summary.txt               # Array job summary
│   ├── *_stats.txt               # Job statistics
│   └── *.out                     # Task outputs
├── startup_benchmark.log         # Startup benchmark log
├── imports_login_benchmark.log   # Import benchmark log (login)
├── imports_compute_benchmark.log # Import benchmark log (compute)
└── array_jobs_benchmark.log      # Array job benchmark log
```

## Next Steps

After running benchmarks:

1. Review `BENCHMARK_SUMMARY.md` in results directory
2. Compare performance metrics against requirements
3. Make go/no-go decision for migration (Task 16.2)
4. If approved, proceed with Phase 7 (Migration Execution)

## Notes

- Benchmarks run 10 iterations each for statistical significance
- Results include mean, standard deviation, min, and max
- Array jobs test both small (10 tasks) and medium (100 tasks) scales
- All benchmarks can be run on login nodes except compute node import tests
- Benchmarks are non-destructive and safe to run multiple times
