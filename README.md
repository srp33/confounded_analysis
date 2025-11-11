# Batch Effect Correction Analysis Pipeline

Analysis pipeline for batch effect correction in gene expression data. Implements the gmm_adjust method alongside deep learning and statistical approaches, with comprehensive evaluation frameworks and flexible environment management.

## Overview

This codebase provides tools for correcting batch effects in gene expression data, with particular focus on the gmm_adjust method and rigorous testing across multiple datasets and evaluation metrics.

**Environment Management:** The pipeline supports two execution approaches:
- **uv/rix (Recommended)**: Fast, native environment management with uv for Python and Nix for R
- **Apptainer (Legacy)**: Container-based execution for HPC clusters

> **Migration Status (Nov 2025)**: Transitioning from Apptainer to uv/rix. Both systems are fully operational during the migration period. See [Migration Timeline](#migration-timeline) below.

## Components

- **Batch Correction**: GMM-based methods (gmm_adjust, gmm_adjust_nonlinear), deep learning approaches (AutoClass, ICVAE, VFAE, Wasserstein), and statistical methods (ComBat, MNN, Seurat, LIGER)
- **Evaluation**: Classification metrics, statistical measures, and biological signal preservation assessment
- **Environment Management**: Native execution with uv/rix or containerization with Docker/Apptainer
- **Data Processing**: Automated pipeline with caching
- **Dataset Support**: 16+ datasets for comprehensive method validation

## System Requirements

### uv/rix (Recommended)
- **Python**: 3.10+ (system Python or module)
- **uv**: Install from [astral.sh/uv](https://docs.astral.sh/uv/)
- **Nix**: Rootless installation via nix-user-chroot (provided in `/grphome/grp_batch_effects/nix/`)
- **Storage**: ~20-30GB in group home for shared environments
- **Cluster**: User namespaces enabled (verified on BYU RC RHEL 9.4)

### Apptainer (Legacy)
- **Apptainer**: 1.0+ (typically available via module system)
- **Storage**: ~41GB for container images
- **Cluster**: Standard HPC environment with SLURM

### Both Approaches
- **SLURM**: For job scheduling (optional, can run interactively)
- **Group Access**: `grp_batch_effects` group membership for shared resources
- **Storage Quotas**: 2 TiB home, 2 TiB group home (BYU RC)

## Quick Start

### uv/rix (Recommended - Native Environments)

**First-time setup:**
```bash
# 1. Initialize environment paths
source environments/init_env.sh

# 2. Install Python environment (one-time, ~30-60 seconds)
cd environments/python
uv sync

# 3. R environment is ready to use (Nix packages are shared)
# No additional setup needed - packages are fetched on first use
```

**Running scripts:**
```bash
# Interactive shell with both Python and R
environments/run_with_env.sh shell

# Run Python scripts
environments/run_with_env.sh scripts/adjust/autoclass.py

# Run R scripts
environments/run_with_env.sh scripts/adjust/gmm_adjust.R

# Submit to SLURM scheduler
environments/run_with_env.sh --sbatch --time 02:00:00 --mem 64G scripts/adjust/gmm_adjust.R

# ComBat-seq environment (Bioconductor 3.11)
environments/run_with_combat_env.sh scripts/evaluations/combat_seq/ComBat_seq.R
```

**Performance:** Environment activation <500ms (vs ~2s for containers), updates in seconds (vs hours for container rebuilds).

### Apptainer (Legacy - Container-based)

```bash
# Interactive shell
./run_in_apptainer.sh shell

# Run specific scripts
./run_in_apptainer.sh scripts/all.sh

# Submit to SLURM scheduler
./run_in_apptainer.sh --sbatch scripts/evaluations/robustifying/code/3_real_data_pipe.R
```

### Docker (Local Development)
```bash
./run_docker.sh
```

## Pipeline Overview

The analysis pipeline consists of four main phases orchestrated by `scripts/all.sh`:

1. **Data Preparation** - Dataset acquisition, processing, and combination generation
2. **Batch Effect Adjustment** - Apply correction methods including gmm_adjust and comparison methods
3. **Evaluation** - Assessment using classification metrics, statistical measures, and biological validation
4. **Visualization** - Generate plots and reports for method comparison

Individual phases can be enabled/disabled for targeted execution. The pipeline includes caching to avoid redundant computations.

See [scripts/README.md](scripts/README.md) for pipeline documentation and configuration options.

## Supported Data

16+ datasets from OSF, Google Drive, and Refinebio including cancer studies (GSE19615, GSE20194, METABRIC) and platform comparisons (GSE96058_HiSeq/NextSeq).

Processing includes download, format standardization, gene annotation, quality control, and pairwise dataset combination generation.

See [scripts/prepdata/README.md](scripts/prepdata/README.md) for data preparation details and [data/README.md](data/README.md) for data structure information.

## Environment Management

### uv/rix Environments (Recommended)

**Architecture:**
- **Python**: Managed by [uv](https://github.com/astral-sh/uv) - ultra-fast package manager with lock files
- **R**: Managed by [rix](https://github.com/b-rodrigues/rix) - Nix-based reproducible R environments
- **Execution**: Native filesystem access (no containers, no bind mounts)
- **Storage**: Shared environments in `/grphome/grp_batch_effects/` for all group members

**Setup:**
```bash
# 1. Source environment initialization
source environments/init_env.sh

# 2. Install Python dependencies (one-time)
cd environments/python
uv sync  # Creates shared .venv in /grphome/grp_batch_effects/environments/python/.venv/

# 3. R environment ready (Nix packages fetched on first use)
# Main analysis: environments/r/batch-effects.nix (Bioconductor 3.21)
# ComBat-seq: environments/r/combatseq.nix (Bioconductor 3.11)
```

**Usage:**
```bash
# Auto-detect script type and activate appropriate environment
environments/run_with_env.sh <script.py|script.R|script.sh>

# SLURM integration
environments/run_with_env.sh --sbatch [sbatch-flags] <script>

# Interactive shell
environments/run_with_env.sh shell

# ComBat-seq environment
environments/run_with_combat_env.sh <script.R>
```

**Benefits:**
- ⚡ Fast startup: <500ms (vs ~2s for containers)
- 🔄 Quick updates: Seconds (vs hours for container rebuilds)
- 📁 Simple paths: Native filesystem access (no bind mounts)
- 💾 Storage efficient: Shared packages across all users
- 🔒 Reproducible: Lock files ensure consistent dependencies

See [environments/README.md](environments/r/README.md) for detailed setup and troubleshooting.

### Apptainer Environments (Legacy)

**Architecture:**
- Three-stage container build (base → fast → annotations)
- Bind mounts for data and scripts
- Group permissions for shared access

**Usage:**
```bash
# Interactive shell
./run_in_apptainer.sh shell

# Run pipeline
./run_in_apptainer.sh scripts/all.sh

# Submit to SLURM scheduler
./run_in_apptainer.sh --sbatch scripts/evaluations/robustifying/code/3_real_data_pipe.R
```

See [apptainer/README.md](apptainer/README.md) for Apptainer documentation including:
- Three-stage build system
- SLURM integration and job submission
- Group permissions setup (`grp_batch_effects`)
- Performance tuning and troubleshooting

### Docker (Local Development)
```bash
# Run container
./run_docker.sh

# Interactive development
docker run -it --rm -v $(pwd):/workspace batch-effects-pipeline bash
```

Single-stage build with Python (Miniforge3), R (Bioconductor 3.21), and dependencies.

## Migration Timeline

**Current Status (November 2025):** Transitioning from Apptainer to uv/rix environments.

### Completed Phases ✅
- **Phase 1 (Investigation):** Nix and uv installations verified on BYU RC cluster
- **Phase 2 (Python Setup):** Python environment production-ready with uv
- **Phase 4 (Init Script):** Environment initialization working (`init_env.sh`)

### In Progress 🔄
- **Phase 3 (R Setup):** R environment created, testing in progress
- **Phase 6 (Documentation):** Updating documentation and migration guides

### Upcoming ⏳
- **Phase 5 (Validation):** Scientific reproducibility validation
- **Phase 7 (Rollout):** Gradual rollout to group members
- **Phase 8 (Production):** Full production deployment

### Progress Summary
- **Investigation:** 80% complete (4/5 tasks)
- **Python Environment:** 100% complete
- **R Environment:** 20% complete (1/5 tasks)
- **Execution Wrappers:** 100% complete
- **Overall:** ~45% complete

**Both systems are fully operational during migration.** Users can choose their preferred execution method:
- Use `environments/run_with_env.sh` for uv/rix (recommended for new work)
- Use `./run_in_apptainer.sh` for Apptainer (stable, fully tested)

See [.kiro/specs/uv-rix-migration/](..kiro/specs/uv-rix-migration/) for detailed migration plan and status.

## Documentation Navigation

## Documentation Structure

### Core Documentation
- [scripts/README.md](scripts/README.md) - Pipeline documentation, execution control, and configuration
- [data/README.md](data/README.md) - Data structure, organization, and management

### Method Documentation  
- [scripts/adjust/README.md](scripts/adjust/README.md) - Batch correction methods (deep learning, statistical, GMM)
- [scripts/evaluations/README.md](scripts/evaluations/README.md) - Evaluation framework and metrics
- [scripts/prepdata/README.md](scripts/prepdata/README.md) - Data preparation and processing pipeline

### Environment Documentation
- [environments/r/README.md](environments/r/README.md) - uv/rix setup, usage, and troubleshooting (recommended)
- [apptainer/README.md](apptainer/README.md) - Apptainer setup, SLURM integration, and HPC deployment (legacy)

### Navigation
- Getting Started: Start here → [scripts/README.md](scripts/README.md) for pipeline details
- Environment Setup: [environments/r/README.md](environments/r/README.md) for uv/rix (recommended) or [apptainer/README.md](apptainer/README.md) for containers
- Method Selection: [scripts/adjust/README.md](scripts/adjust/README.md) for batch correction methods
- Results Analysis: [scripts/evaluations/README.md](scripts/evaluations/README.md) for evaluation metrics
- Data Management: [scripts/prepdata/README.md](scripts/prepdata/README.md) + [data/README.md](data/README.md)
- HPC Deployment: [environments/r/README.md](environments/r/README.md) for native execution or [apptainer/README.md](apptainer/README.md) for containers

## Available Methods

### Batch Correction Methods
- **Gaussian Mixture Models**: gmm_adjust (primary method), gmm_adjust_nonlinear, gmm_global_simple
- **Deep Learning**: AutoClass, ICVAE, VFAE, Wasserstein adversarial approaches
- **Statistical**: ComBat, MNN, Seurat integration, LIGER, limma, quantile normalization

### Evaluation Metrics
- **Classification**: Batch vs. biological signal separation, ER status prediction
- **Statistical**: MMD, MSE, mutual information, feature importance analysis
- **Quality Control**: BatchQC, pathway analysis, dimensionality reduction

See [scripts/adjust/README.md](scripts/adjust/README.md) for method descriptions and [scripts/evaluations/README.md](scripts/evaluations/README.md) for evaluation framework details.

### Environment Comparison

| Feature | uv/rix (Recommended) | Apptainer (Legacy) |
|---------|---------------------|-------------------|
| **Startup Time** | <500ms | ~2s |
| **Environment Updates** | Seconds | Hours (rebuild) |
| **Path Handling** | Native filesystem | Bind mounts |
| **Storage** | ~18-27GB shared | ~41GB images |
| **Setup Time** | <5 minutes | 2+ hours |
| **Multi-user** | Shared packages | Shared images |
| **Flexibility** | Per-script environments | Single environment |
| **Status** | Production-ready | Stable, fully tested |

**When to use uv/rix:**
- New projects and analyses
- Rapid prototyping and development
- Frequent environment updates
- Need for fast job startup

**When to use Apptainer:**
- Existing validated workflows
- Maximum stability required
- Prefer containerization
- During migration validation

### Diagnostic Commands

**uv/rix:**
```bash
# Check environment status
source environments/init_env.sh

# Test Python environment
source /grphome/grp_batch_effects/environments/python/.venv/bin/activate
python -c "import numpy, torch, sklearn; print('Python OK')"

# Test R environment
/grphome/grp_batch_effects/nix/nix-env nix-shell environments/r/batch-effects.nix --run "Rscript -e 'library(tidyverse); print(\"R OK\")'"

# Interactive shell
environments/run_with_env.sh shell

# Check storage usage
du -sh /grphome/grp_batch_effects/environments
du -sh /grphome/grp_batch_effects/.uv_cache
```

**Apptainer:**
```bash
# Check pipeline logs
tail -f outputs/prepdata.log

# Interactive container access
./run_in_apptainer.sh shell
```

See [environments/r/README.md](environments/r/README.md) for uv/rix troubleshooting or [apptainer/README.md](apptainer/README.md) for Apptainer optimization strategies.

## Development and Contributing

### Development Workflow

1. **Method Development**: Add new batch correction methods in `scripts/adjust/`
2. **Evaluation Extensions**: Extend metrics and validation in `scripts/evaluations/`
3. **Data Integration**: Add new data sources in `scripts/prepdata/`
4. **Container Updates**: Modify container definitions in `apptainer/`

### Contributing Guidelines
See method-specific documentation in:
- [scripts/adjust/README.md](scripts/adjust/README.md) - Adding new batch correction methods
- [scripts/evaluations/README.md](scripts/evaluations/README.md) - Extending evaluation frameworks
- [scripts/prepdata/README.md](scripts/prepdata/README.md) - Integrating new data sources
- [environments/r/README.md](environments/r/README.md) - Environment management with uv/rix
- [apptainer/README.md](apptainer/README.md) - Container development (legacy)

### Support and Contact
- Issues: [GitHub Issues](https://github.com/srp33/confounded_analysis/issues) for bug reports and feature requests
- Research Contact: [Piccolo Lab](https://biology.byu.edu/piccolo-lab/contact) for research collaboration

