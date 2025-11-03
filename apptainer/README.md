# Apptainer Container Documentation

> **Navigation**: [← Main README](../README.md) | [Pipeline Documentation →](../scripts/README.md) | [Data Structure →](../data/README.md)

Documentation for Apptainer/Singularity container setup, HPC integration, SLURM usage, and troubleshooting for the batch effect correction analysis pipeline.

## Container Architecture

### Three-Stage Build System

The build system uses a three-stage approach for development efficiency, resource utilization, and maintenance:

#### Stage 1: Base Image (`apptainer_base.def`)
- **Purpose**: Stable foundation with system dependencies and core R packages
- **Contents**: 
  - Ubuntu 22.04 LTS base system
  - R 4.4+ with Bioconductor 3.21
  - Core system libraries and build tools
  - Essential R packages for genomics analysis
  - ccache compilation optimization
- **Rebuild Frequency**: Only when system dependencies change

#### Stage 2: Fast Image (`apptainer_fast.def`)
- **Purpose**: Python environment and specialized scientific packages
- **Contents**:
  - Miniforge3 Python environment
  - Deep learning frameworks (TensorFlow, PyTorch)
  - Scientific computing stack (NumPy, SciPy, scikit-learn)
  - Specialized packages (AIF360, UMAP, openTSNE)
  - AutoClass deep learning framework
- **Rebuild Frequency**: When Python packages or methods are updated

#### Stage 3: Annotations Image (`apptainer_annotations.def`)
- **Purpose**: Additional annotation tools and specialized packages
- **Contents**:
  - Custom annotation packages from external sources
  - Specialized genomics tools
  - Additional R packages for pathway analysis
  - Extended visualization libraries
- **Rebuild Frequency**: When annotation tools are added or updated

### Build Dependencies and Optimization

**Dependency Chain**:
```
Ubuntu 22.04 → Base Image → Fast Image → Annotations Image
     ↓              ↓            ↓              ↓
System Deps    Core R Pkgs   Python Env   Annotations
```

### Complete Build Process

#### Initial Setup and Prerequisites

**System Requirements**:
- **Permissions**: Access to `grp_batch_effects` group and SLURM scheduler

**Environment Setup**:
```bash
# Load required modules
module load apptainer

# Set up build directories
mkdir -p ~/confounded_analysis/grp_batch_effects
mkdir -p ~/confounded_analysis/build_logs

# Initialize environment variables
source ~/confounded_analysis/init_apptainer.sh
```

#### Automated Build Chain

**Complete Automated Build**:
```bash
# Submit entire build chain with dependencies
sbatch build_base_image.sh
# build_base_image.sh automatically submits build_fast_image.sh upon completion
# build_fast_image.sh automatically submits build_annotations_image.sh upon completion

# Monitor all builds
watch 'squeue -u $USER | grep build'
```


### Container Definitions

**Base Container** (`apptainer_base.def`):
- **Base OS**: Ubuntu 22.04 LTS
- **R Environment**: Bioconductor 3.21 with core genomics packages
- **System Dependencies**: Build tools, libraries, and system packages
- **Purpose**: Stable foundation for incremental builds

**Fast Container** (`apptainer_fast.def`):
- **Python Environment**: Miniforge3 with scientific computing stack
- **Deep Learning**: TensorFlow ≥2.0, PyTorch with GPU support
- **Specialized Libraries**: AIF360, UMAP, openTSNE, fairness tools
- **Purpose**: Rapid development and method iteration

**Annotations Container** (`apptainer_annotations.def`):
- **Additional Tools**: Specialized annotation and analysis tools
- **Extended Libraries**: Additional R and Python packages
- **Purpose**: Analysis capabilities

## Usage Examples

### Basic Usage

**Interactive Shell**:
```bash
# Start interactive shell
./run_in_apptainer.sh shell

# Interactive shell with custom image
./run_in_apptainer.sh --image-path /path/to/custom.sif shell
```

### SLURM Integration

**Basic SLURM Submission**:
```bash
# Submit script to SLURM scheduler
./run_in_apptainer.sh --sbatch scripts/evaluations/robustifying/code/3_real_data_pipe.R

# Submit with custom resources
./run_in_apptainer.sh --sbatch --time 02:00:00 --mem 64G scripts/adjust/autoclass.py
```

**SLURM Options**:
```bash
# Submit with specific partition and QOS
./run_in_apptainer.sh --sbatch --partition gpu --qos gpu --gres gpu:1 scripts/adjust/icvae.py

# Submit array jobs
./run_in_apptainer.sh --sbatch --array 1-10 scripts/adjust/gmm_adjust.R
```

### Example SLURM Workflows

**Complete ComBat-seq Workflow**:
```bash
# Submit complete ComBat-seq workflow
sbatch grp_batch_effects/slurm_scripts/combat_seq/combat_seq_complete_workflow.sh

# Monitor job progress
squeue -u $USER | grep combat_seq
```

**Robustification Testing**:
```bash
# Run robustification testing
sbatch grp_batch_effects/slurm_scripts/robustifying/robustifying_complete_workflow.sh

# Check job dependencies
squeue -u $USER --format="%.18i %.9P %.50j %.8u %.8T %.10M %.9l %.6D %R"
```

### Shared Resource Management

**Container Images**:
- Location: `~/groups/grp_batch_effects/remove-batch-effects*.sif`
- Access control: Read access for all group members, write access for administrators
- Version management: Centralized updates with backward compatibility
- Backup strategy: Regular snapshots and version archiving

**Data Directories**:
```bash
# Shared data structure
~/groups/grp_batch_effects/
├── data/                    # Shared datasets and cache
│   ├── .cache/             # Shared computation cache
│   └── datasets/           # Processed datasets
├── outputs/                # Shared analysis outputs
│   ├── figures/           # Generated plots and visualizations
│   ├── metrics/           # Evaluation results
│   └── tables/            # Summary statistics
├── remove-batch-effects-base.sif
├── remove-batch-effects-fast.sif
├── remove-batch-effects.sif
└── slurm_outputs/         # SLURM job logs and outputs
```

**Access Patterns**:
```bash
# Individual user workspace (private)
~/confounded_analysis/

# Shared group resources (collaborative)
~/groups/grp_batch_effects/

# Bind mount configuration (automatic via init_apptainer.sh)
export APPTAINER_BINDPATH="$SHARED_DIR/data:/data,$ANALYSIS_DIR/apptainer:/apptainer,$SHARED_DIR/outputs:/outputs,$SCRIPTS_DIR:/scripts,$SHARED_DIR:$ANALYSIS_DIR/grp_batch_effects_folders_in_apptainer_are_located_at_root"
```


**Next step: Integration with Workflow Managers**:
```bash
# 1. Nextflow integration
nextflow run pipeline.nf \
    --container_image $APPTAINER_IMAGE \
    --executor slurm \
    --queue compute

# 2. Snakemake integration
snakemake --use-singularity \
    --singularity-args "--bind /data:/data" \
    --cluster "sbatch --mem={resources.mem_mb} --time={resources.time}"

# 3. CWL integration
cwltool --singularity --tmp-outdir-prefix /tmp/ workflow.cwl inputs.yml
```

## Troubleshooting Guide

### Build-Time Issues

**Memory-Related Build Failures**:
```bash
# Symptom: Build killed with exit code 137 (SIGKILL)
# Cause: Out of memory during compilation

# Solution 1: Increase SLURM memory allocation
sbatch --mem=64G build_base_image.sh

# Solution 2: Use build nodes with more memory
sbatch --partition=highmem --mem=128G build_base_image.sh

# Solution 3: Enable swap if available
swapon --show  # Check current swap
```

**Network and Dependency Issues**:
```bash
# Symptom: Package download failures, DNS resolution errors
# Diagnosis commands:
getent hosts cloud.r-project.org
wget -q --spider https://conda-forge.org
nslookup cran.r-project.org

# Solutions:
# 1. Check network connectivity
ping 8.8.8.8

# 2. Verify DNS configuration
cat /etc/resolv.conf

# 3. Use alternative package sources
# Edit apptainer_base.def to use different CRAN mirror
```

**Permission and File System Issues**:
```bash
# Symptom: Permission denied during build
# Check build directory permissions
ls -la ~/confounded_analysis/grp_batch_effects/

# Fix permissions
chmod 755 ~/confounded_analysis/grp_batch_effects/
mkdir -p /var/tmp/$USER/apptainer_build_tmp
chmod 755 /var/tmp/$USER/apptainer_build_tmp

# Verify fakeroot capability
apptainer build --help | grep fakeroot
```

## Configuration

### Custom Container Builds

**Custom Definition Files**:
```bash
# Create custom container definition
cp apptainer_fast.def my_custom.def
# Edit my_custom.def with additional packages

# Build custom container
apptainer build --fakeroot my_custom.sif my_custom.def
```

---

> **Navigation**: [← Main README](../README.md) | [Pipeline Documentation →](../scripts/README.md) | [Data Structure →](../data/README.md)