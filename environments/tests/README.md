# Execution Wrapper Tests (Task 10.4)

This directory contains comprehensive tests for the `run_with_env.sh` execution wrapper and related scripts.

## Test Files

### Test Scripts

1. **test_python_script.py** - Tests Python environment execution
   - Verifies Python version and executable
   - Tests package imports (NumPy, Pandas, scikit-learn)
   - Tests command-line argument handling
   - Performs simple computations

2. **test_r_script.R** - Tests R environment execution
   - Verifies R version and configuration
   - Tests package loading (dplyr, ggplot2, data.table)
   - Tests command-line argument handling
   - Performs data manipulations

3. **test_mixed_env.R** - Tests mixed Python/R environment
   - Tests R with Python integration via reticulate
   - Verifies RETICULATE_PYTHON environment variable
   - Tests Python package imports from R
   - Tests R data manipulation with dplyr

4. **test_shell_script.sh** - Tests shell script execution
   - Verifies both Python and R are available
   - Tests environment variables
   - Tests package availability in both languages

### Test Runner

**run_all_tests.sh** - Comprehensive test suite that runs all tests

## Running Tests

### Run All Tests

```bash
cd environments/tests
./run_all_tests.sh
```

### Run Individual Tests

#### Python Script Test
```bash
cd environments
./run_with_env.sh tests/test_python_script.py
./run_with_env.sh tests/test_python_script.py arg1 arg2
```

#### R Script Test
```bash
cd environments
./run_with_env.sh tests/test_r_script.R
./run_with_env.sh --r-env combatseq tests/test_r_script.R
```

#### Mixed Environment Test
```bash
cd environments
./run_with_env.sh --full-env tests/test_mixed_env.R
```

#### Shell Script Test
```bash
cd environments
./run_with_env.sh tests/test_shell_script.sh
```

### Test Interactive Shell
```bash
cd environments
./run_with_env.sh shell

# Inside the shell, test:
python --version
R --version
python -c "import numpy; print(numpy.__version__)"
Rscript -e "library(dplyr); cat('dplyr loaded\n')"
exit
```

### Test SLURM Job Submission

**Note:** These tests require access to a SLURM cluster.

```bash
cd environments

# Python job
./run_with_env.sh --sbatch --time 00:10:00 --mem 4G tests/test_python_script.py

# R job
./run_with_env.sh --sbatch --time 00:10:00 --mem 4G tests/test_r_script.R

# Mixed environment job
./run_with_env.sh --sbatch --time 00:10:00 --mem 4G --full-env tests/test_mixed_env.R

# With custom resources
./run_with_env.sh --sbatch --time 02:00:00 --mem 64G --cpus-per-task 8 tests/test_python_script.py

# With array jobs
./run_with_env.sh --sbatch --array 1-10 --time 00:05:00 --mem 2G tests/test_python_script.py
```

## Test Coverage

### Task 10.4 Sub-tasks

- ✓ **Test Python script execution** - `test_python_script.py`
- ✓ **Test R script execution** - `test_r_script.R`
- ✓ **Test mixed environment execution** - `test_mixed_env.R`, `test_shell_script.sh`
- ✓ **Test SLURM job submission** - `run_all_tests.sh` (Category 6)
- ✓ **Test with various resource specifications** - `run_all_tests.sh` (Category 7)

### Additional Test Categories

1. **Python Script Execution** (Category 1)
   - Basic execution
   - With command-line arguments

2. **R Script Execution** (Category 2)
   - Basic execution
   - With command-line arguments
   - With ComBat-seq environment

3. **Mixed Environment Execution** (Category 3)
   - Full environment flag
   - Shell scripts with both environments
   - Reticulate integration

4. **Direct Command Execution** (Category 4)
   - Python commands
   - R commands

5. **Interactive Shell** (Category 5)
   - Manual verification required

6. **SLURM Job Submission** (Category 6)
   - Script generation
   - Job submission (requires cluster)

7. **Resource Specifications** (Category 7)
   - Various time/memory combinations
   - CPU and node specifications
   - Partition selection

8. **Error Handling** (Category 8)
   - Missing files
   - Invalid environments

## Expected Results

All tests should pass with the following indicators:

- ✓ Green checkmarks for passed tests
- Detailed output showing:
  - Environment variables set correctly
  - Packages loaded successfully
  - Computations executed correctly
  - Command-line arguments passed through

## Troubleshooting

### Python Environment Not Found

```bash
cd environments/python
uv sync
```

### R Environment Not Found

Check that the Nix environment files exist:
```bash
ls -l environments/r/batch-effects.nix
ls -l environments/r/combatseq.nix
```

### Nix Wrapper Not Found

Verify the Nix installation:
```bash
ls -l /grphome/grp_batch_effects/nix/nix-env
```

### Package Import Failures

For Python:
```bash
source /grphome/grp_batch_effects/environments/python/.venv/bin/activate
python -c "import numpy, pandas, sklearn"
```

For R:
```bash
/grphome/grp_batch_effects/nix/nix-env nix-shell environments/r/batch-effects.nix --run "Rscript -e 'library(dplyr)'"
```

## Requirements Verification

This test suite verifies **Requirement 4.5** from the requirements document:

> THE Execution_Wrapper SHALL maintain backward compatibility with existing script invocation patterns during migration

The tests ensure:
- Scripts can be executed with simple commands
- Command-line arguments are passed through correctly
- Environment variables are set properly
- Both Python and R environments work correctly
- SLURM integration functions as expected
- Various resource specifications are handled correctly
