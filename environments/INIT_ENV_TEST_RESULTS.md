# init_env.sh Test Results

## Summary

The `init_env.sh` script has been successfully implemented and tested on both login nodes and SLURM compute nodes. All tests passed successfully.

## Test Date

November 7, 2025

## Tests Performed

### 1. Login Node Tests (test_init_env.sh)

**Status:** ✅ PASSED

**Test Coverage:**
- ✅ All required environment variables are set correctly
- ✅ All paths use explicit notation (`/grphome/` or `${HOME}/`)
- ✅ All required directories are created
- ✅ Group permissions are configured correctly (775 with setgid)
- ✅ Python environment is added to PATH
- ✅ `check_storage_quota()` function works correctly
- ✅ `verify_environment()` function works correctly
- ✅ Verbose mode (`VERBOSE=1`) works correctly

**Key Results:**
- Storage usage: 596M in user home, 758G in group home (well within 2 TiB quotas)
- Python environment: Python 3.12.2 with NumPy 2.3.4
- All critical components verified

### 2. SLURM Compute Node Tests (test_init_env_slurm.sh)

**Status:** ✅ PASSED

**Job Details:**
- Job ID: 8367979
- Node: m8-17-11
- Exit Code: 0
- Elapsed Time: 2 seconds

**Test Coverage:**
- ✅ Script sources successfully on compute nodes
- ✅ All environment variables are accessible
- ✅ All directories are accessible from compute nodes
- ✅ Python environment works (Python 3.12.2, NumPy 2.3.4)
- ✅ Nix wrapper is accessible and executable
- ✅ Write access to shared directories confirmed
- ✅ Functions are defined and available

## Environment Variables Verified

All required environment variables are set correctly:

| Variable | Value | Status |
|----------|-------|--------|
| SHARED_DIR | /grphome/grp_batch_effects | ✅ |
| ANALYSIS_DIR | ${HOME}/confounded_analysis | ✅ |
| SCRIPTS_DIR | ${HOME}/confounded_analysis/scripts | ✅ |
| DATA_DIR | /grphome/grp_batch_effects/data | ✅ |
| OUTPUTS_DIR | /grphome/grp_batch_effects/outputs | ✅ |
| SCRATCH_DIR | /grphome/grp_batch_effects/nobackup/autodelete | ✅ |
| ARCHIVE_DIR | /grphome/grp_batch_effects/nobackup/archive | ✅ |
| PYTHON_ENV_SPEC | ${HOME}/confounded_analysis/environments/python | ✅ |
| PYTHON_ENV | /grphome/grp_batch_effects/environments/python/.venv | ✅ |
| UV_CACHE_DIR | /grphome/grp_batch_effects/.uv_cache | ✅ |
| UV_PROJECT_ENVIRONMENT | /grphome/grp_batch_effects/environments/python/.venv | ✅ |
| PYTHON | /grphome/grp_batch_effects/environments/python/.venv/bin/python | ✅ |
| PYTHON3 | /grphome/grp_batch_effects/environments/python/.venv/bin/python | ✅ |
| RETICULATE_PYTHON | /grphome/grp_batch_effects/environments/python/.venv/bin/python | ✅ |
| R_ENV_BATCH_EFFECTS | ${HOME}/confounded_analysis/environments/r/batch-effects.nix | ✅ |
| R_ENV_COMBATSEQ | ${HOME}/confounded_analysis/environments/r/combatseq.nix | ✅ |
| R_LIBS_USER | /grphome/grp_batch_effects/environments/r/r-libs | ✅ |
| NIX_WRAPPER | /grphome/grp_batch_effects/nix/nix-env | ✅ |
| NIX_PATH | nixpkgs=channel:nixos-23.11 | ✅ |

## Directories Created

All required directories were successfully created:

- ✅ /grphome/grp_batch_effects/data
- ✅ /grphome/grp_batch_effects/data/.cache
- ✅ /grphome/grp_batch_effects/outputs/figures
- ✅ /grphome/grp_batch_effects/outputs/metrics
- ✅ /grphome/grp_batch_effects/outputs/tables
- ✅ /grphome/grp_batch_effects/environments/python
- ✅ /grphome/grp_batch_effects/environments/r/r-libs
- ✅ /grphome/grp_batch_effects/.uv_cache
- ✅ /grphome/grp_batch_effects/nobackup/autodelete (scratch)
- ✅ /grphome/grp_batch_effects/nobackup/archive

## Permissions

Group permissions are correctly configured:
- Directories: 775 (rwxrwxr-x)
- Files: 664 (rw-rw-r--)
- Setgid bit set on shared directories for group ownership inheritance

## Functions Implemented

### check_storage_quota()

Displays storage usage across all BYU RC storage tiers:
- User Home (2 TiB quota, backed up)
- Group Home (2 TiB quota, backed up)
- Group Scratch (20 TiB quota, 12-week auto-delete)
- Group Archive (20 TiB quota, long-term storage)
- Warns if approaching quota limits (>90% usage)

**Usage:**
```bash
source environments/init_env.sh
check_storage_quota
```

Or:
```bash
SHOW_STORAGE=1 source environments/init_env.sh
```

### verify_environment()

Verifies that all critical components are present:
- Python environment and version
- R environment specifications
- Nix wrapper
- Data and output directories

**Usage:**
```bash
source environments/init_env.sh
verify_environment
```

Or:
```bash
VERIFY=1 source environments/init_env.sh
```

## Usage Examples

### Basic Usage
```bash
source environments/init_env.sh
```

### Verbose Mode
```bash
VERBOSE=1 source environments/init_env.sh
```

### Show Storage Usage
```bash
SHOW_STORAGE=1 source environments/init_env.sh
```

### Verify Environment
```bash
VERIFY=1 source environments/init_env.sh
```

### Combined Options
```bash
VERBOSE=1 SHOW_STORAGE=1 VERIFY=1 source environments/init_env.sh
```

### In SLURM Jobs
```bash
#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --mem=32G

# Initialize environment
source ${HOME}/confounded_analysis/environments/init_env.sh

# Run your script
python scripts/my_analysis.py
```

## Requirements Met

This implementation satisfies all requirements from task 9.1:

- ✅ Set all environment variables (paths, caches)
- ✅ Use explicit paths (`/grphome/`, `${HOME}/`)
- ✅ Create necessary directories
- ✅ Set group permissions (775 with setgid)
- ✅ Add storage quota checking function
- ✅ Requirements: 4.4, 5.4, 5.5

## Next Steps

Task 9 (Create environment initialization script) is now complete. The script is ready for production use.

Next recommended tasks:
- Task 10.1: Create `run_with_env.sh` execution wrapper
- Task 8.3: Build and test batch-effects R environment
- Task 4: Benchmark performance vs Apptainer

## Notes

- The storage quota checking function can be slow on large directories (uses `du` commands)
- For SLURM jobs, consider skipping storage checks to reduce startup time
- All paths use explicit notation as required (no tilde expansion)
- Group permissions ensure all members can access shared environments
- Setgid bit ensures new files inherit group ownership
