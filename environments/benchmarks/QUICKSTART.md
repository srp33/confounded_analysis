# Benchmark Quick Start Guide

## Run All Benchmarks (Recommended)

```bash
cd ~/confounded_analysis/environments/benchmarks
bash run_all_benchmarks.sh
```

This will run all benchmarks sequentially and generate a comprehensive summary report.

**Time required:** ~30-60 minutes (includes waiting for SLURM jobs)

## Run Individual Benchmarks

### 1. Startup Times (~5 minutes)

```bash
bash benchmark_startup.sh
```

Results: `startup_times.csv`

### 2. Import Times - Login Node (~10 minutes)

```bash
bash benchmark_imports.sh login
```

Results: `import_times.csv`

### 3. Import Times - Compute Node (~15 minutes)

```bash
sbatch --time=01:00:00 --mem=4G --wrap="bash $(pwd)/benchmark_imports.sh compute"
```

Results: Appended to `import_times.csv`

### 4. Array Jobs (~20-30 minutes)

```bash
bash benchmark_array_jobs.sh
```

Results: `array_job_results/` directory

## View Results

```bash
# View startup times
cat startup_times.csv

# View import times
cat import_times.csv

# View array job summary
cat array_job_results/summary.txt

# View comprehensive summary (after running all benchmarks)
cat results_*/BENCHMARK_SUMMARY.md
```

## Prerequisites

- ✅ Python environment set up (Task 7)
- ⚠️ R environment built (Task 8.3) - Optional, R benchmarks will be skipped if not available
- ✅ Apptainer image available
- ✅ Access to SLURM

## Expected Output

Each benchmark will:
1. Run 10 iterations for statistical significance
2. Display progress in real-time
3. Calculate mean, standard deviation, min, and max
4. Save results to CSV files
5. Generate summary reports

## Troubleshooting

**"R environment not found"**
- R benchmarks will be skipped automatically
- To enable R benchmarks, complete Task 8.3 first

**"Python environment not found"**
- Run: `cd environments/python && uv sync`

**SLURM jobs pending**
- Check queue: `squeue -u $USER`
- Jobs may take time to start depending on cluster load

## Next Steps

After benchmarks complete:
1. Review results in `results_*/BENCHMARK_SUMMARY.md`
2. Compare against requirements (job startup < 500ms)
3. Make go/no-go decision for migration
