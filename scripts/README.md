# Pipeline Documentation

> **Navigation**: [← Main README](../README.md) | [Batch Correction Methods →](adjust/README.md) | [Evaluation Framework →](evaluations/README.md) | [Data Preparation →](prepdata/README.md)

## Pipeline Architecture

The analysis pipeline consists of four main phases.

### 1. Data Preparation Phase (`scripts/prepdata/all.sh`)

**Purpose**: Dataset acquisition, processing, and preparation for batch correction analysis.

**Components**:
- Multi-source downloads from OSF, Google Drive, and Refinebio
  - OSF datasets: GSE19615, GSE20194, GSE20271, GSE23720, GSE25055, GSE25065, GSE31448, GSE45255, GSE58644, GSE62944_Tumor, GSE76275, GSE81538, GSE96058_HiSeq, GSE96058_NextSeq, METABRIC
  - Google Drive datasets: GSE115577, GSE123845, GSE163882
- Dataset processing: Format standardization, gene ID mapping

**Output Structure**:
- `/data/gold/[dataset]/unadjusted.csv` - Processed individual datasets
- `/data/paired_datasets/` - Combined dataset pairs for cross-study analysis
- `/outputs/prepdata.log`, `/outputs/prepdata2.log` - Processing logs

See [prepdata/README.md](prepdata/README.md) for details.

### 2. Batch Effect Adjustment Phase (`scripts/adjust/all.sh`)

**Purpose**: Apply various batch correction methods to prepared datasets with optimized resource allocation.

**Method Categories**:
- **Deep Learning Methods**: AutoClass, ICVAE, VFAE, Wasserstein adversarial approaches
- **Statistical Methods**: ComBat, quantile normalization, MNN, Seurat integration, LIGER, limma
- **Gaussian Mixture Models**: Multiple GMM variants with nonlinear extensions
- **Ranking Methods**: ranked1, ranked2, ranked_batch (currently active)

**Execution Strategy**:
- **Parallel Methods**: Lightweight methods run simultaneously across datasets
- **Sequential Methods**: Memory-intensive methods run one dataset at a time
- **Target-Preserving Methods**: ComBat, FairAdapt, Limma with biological signal preservation

See [adjust/README.md](adjust/README.md) for details.

### 3. Evaluation Phase (`scripts/evaluations/all.sh`)

**Purpose**: Multi-dimensional assessment of batch correction effectiveness using diverse metrics.

**Evaluation Categories**:
- **Classification Assessment**: Batch vs. biological signal separation, ER status prediction
- **Statistical Metrics**: MMD, MSE, mutual information, feature importance analysis
- **Specialized Analyses**: BatchQC quality control, ESR1 pathway analysis, dimensionality reduction
- **Robustness Testing**: Ensemble learning validation, cross-modal performance assessment

See [evaluations/README.md](evaluations/README.md) for details.

### 4. Visualization and Reporting Phase

**Purpose**: Generate visualizations and reports for method comparison and results interpretation.

**Components**:
- **Performance Plots**: Method comparison visualizations
- **Dimensionality Reduction**: t-SNE, UMAP, PCA plots
- **Statistical Reports**: LaTeX-formatted tables and analysis reports
- **Classification Results**: Confusion matrices and performance summaries

## Pipeline Execution Control

### Master Pipeline Script (`scripts/all.sh`)

The main orchestration script controls the entire pipeline execution.

**Key Features**:
- **Error Handling**: `set -e` ensures pipeline stops on first error
- **Environment Setup**: PYTHONPATH configuration for module imports
- **Logging**: Phase-specific log redirection to `/outputs/` directory
- **Modular Control**: Individual phases can be enabled/disabled via commenting

### Execution Modes


**Container-Based Execution**:
```bash
# Docker
./run_docker.sh

# Interactive container access
./run_in_apptainer.sh shell

# Direct pipeline execution
./run_in_apptainer.sh scripts/all.sh

# SLURM job submission
./run_in_apptainer.sh --sbatch --time 04:00:00 --mem 128G scripts/adjust/adjustR_data.sh
```

## Configuration Options

### Resource Allocation and Optimization

**Memory Management**:
```bash
# Increase stack size for containerized environments
ulimit -s unlimited

# SLURM resource allocation
./run_in_apptainer.sh --sbatch --time 04:00:00 --mem 128G --cpus-per-task 16
```

## Caching System

### Cache Architecture

The pipeline implements a caching system to avoid redundant computations and enable incremental processing.

**Cache Components**:
- **HashCache Class**: Content-based validation using MD5 hashes
- **DataFrameCache Class**: In-memory caching for frequently accessed data
- **File-Based Locking**: Concurrent access protection for multi-process execution

**Cache Locations**:
```bash
data/.cache/                    # Main cache directory
├── gmm_cache/                 # GMM method-specific cache
├── gdown/                     # Download cache
├── R/                         # R package cache
├── classify_hashes/           # Classification cache
└── [method]_cache.json        # Method-specific hash files
```

**Cache Clearing Strategies**:
```bash
# Clear all caches (force complete rerun)
rm -rf data/.cache/

# Clear method-specific caches
rm -rf data/.cache/gmm_cache/
rm -rf data/.cache/R/

# Clear download caches only
rm -rf data/.cache/gdown/
```

**Hash-Based Validation**:
- Input file changes automatically invalidate cache
- Content-based hashing prevents false cache hits
- Incremental processing for partial result completion

### Performance Optimization Features

**Incremental Processing**:
- Automatic detection of completed vs. incomplete runs
- "Top-up" functionality for partial result sets

## Performance Monitoring and Logging

### Logging Infrastructure

**Phase-Specific Logs**:
```bash
/outputs/prepdata.log          # Data preparation phase
/outputs/prepdata2.log         # Combination generation
/outputs/adjust.log            # Batch correction phase
/outputs/metrics.log           # Evaluation metrics
/outputs/figures.log           # Visualization generation
```

**Method-Specific Logs**:
```bash
/outputs/esr1_analysis.log
/outputs/classify.log
/outputs/hist_gradient_er.log
/outputs/feature_importance.log
```

### Method-Level Customization

**Data Preparation Customization** (`scripts/prepdata/all.sh`):
```bash
# Enable/disable download sources
# OSF datasets (comment to disable)
python3 /scripts/prepdata/download_datasets.py --source osf --datasets "$osf_datasets"

# Google Drive datasets (comment to disable)  
# python3 /scripts/prepdata/download_datasets.py --source gdrive --datasets "$gdrive_datasets"

# Control combination generation
python3 /scripts/prepdata/generate_all_combinations.py \
    --csv-files unadjusted.csv \
    --debug \
    --parallel 10  # Adjust parallelism based on resources
```

**Adjustment Method Selection** (`scripts/adjust/adjustR_data.sh`):
```bash
# Lightweight methods (run in parallel)
ADJUSTERS_PARALLEL=(
    "quantile"          # Enable quantile normalization
    # "min_mean"        # Disable min-mean normalization
    "combat"            # Enable ComBat
    "ranked1"           # Enable ranked method 1
    "ranked2"           # Enable ranked method 2
    "ranked_batch"      # Enable batch-aware ranking
)

# Memory-intensive methods (run sequentially)
ADJUSTERS_SEQUENTIAL=(
    # "mnn"             # Disable MNN (memory intensive)
    # "liger"           # Disable LIGER (memory intensive)
)

# Target-preserving methods
ADJUSTERS_TARGET=(
    "combat"            # ComBat with target preservation
    # "fairadapt"       # Disable FairAdapt
    # "limma"           # Disable limma
)
```

**Evaluation Method Selection** (`scripts/evaluations/all.sh`):
```bash
# Enable specific evaluation components
# python /scripts/evaluations/esr1/esr1_analysis.py &> /outputs/esr1_analysis.log

# Classification evaluations
bash /scripts/evaluations/classify_er_mixed_datasets/classify.sh &> /outputs/hist_gradient_er.log

# Statistical evaluations (currently disabled)
# bash /scripts/evaluations/small_evals/mutual_info.sh &> /outputs/mutual_info.log
# bash /scripts/evaluations/small_evals/mse.sh &> /outputs/mse.log
# bash /scripts/evaluations/small_evals/mmd.sh &> /outputs/mmd.log
```

### Dataset Configuration

**Dataset Selection** (`scripts/adjust/adjustR_data.sh`):
```bash
# Active datasets
DATASETS=(
    "gse49711"          # Neuroblastoma study
    "gse20194"          # Breast cancer study  
    "gse24080"          # Pediatric cancer study
    # "gse62944"        # Additional dataset (disabled)
)

# Synthetic datasets for testing
SYNTHETIC_DATASETS=(
    # "2_dims_no_bio_no_batch"
    # "400_dims_yes_bio_yes_batch"
    # "1000_dims_yes_bio_no_batch"
    # "structured_synthetic"
)
```

**Batch Column Mapping**:
```bash
# Configure batch columns for each dataset
declare -A BATCH_COLS
BATCH_COLS["gse49711"]="meta_Sex"              # Sex as batch variable
BATCH_COLS["gse20194"]="meta_Dataset_ID"       # Dataset ID as batch
BATCH_COLS["gse24080"]="meta_batch"            # Explicit batch column
BATCH_COLS["custom_dataset"]="meta_platform"   # Platform as batch variable
```

**Target Preservation Configuration**:
```bash
# Biological signals to preserve during adjustment
declare -A TARGET_COLS
TARGET_COLS["gse49711"]="meta_INSS_Stage_Split_3_4"  # Cancer stage
TARGET_COLS["gse20194"]="meta_er_status"             # ER status
TARGET_COLS["gse24080"]="meta_efs_outcome_label"     # Survival outcome
```

## Pipeline Configuration

### Resource Management Strategies


**SLURM Integration**:
```bash
# Job submission with resource specifications
./run_in_apptainer.sh --sbatch \
    --job-name="batch_correction" \
    --time=08:00:00 \
    --mem=128G \
    --cpus-per-task=16 \
    --partition=gpu \
    scripts/adjust/all.sh

# Array job submission for multiple datasets
./run_in_apptainer.sh --sbatch \
    --array=1-3 \
    --job-name="adjust_array" \
    scripts/adjust/run_single_adjust.sh
```

### Pipeline Customization Patterns

**Custom Method Integration**:
```bash
# Add new adjustment method to adjustR_data.sh
ADJUSTERS_PARALLEL+=(
    "custom_method"     # Add your custom method
)

# Configure method-specific parameters
declare -A METHOD_PARAMS
METHOD_PARAMS["custom_method"]="--param1 value1 --param2 value2"

# Modify run_adjust function to use parameters
run_adjust() {
    local adjuster=$1
    local dataset=$2
    local params="${METHOD_PARAMS[$adjuster]:-""}"
    
    Rscript "$ADJUST_SCRIPT" "$input_file" "$output_file" \
        -a "$adjuster" -b "$batch_col" $params $c_args --debug
}
```

**Custom Evaluation Metrics**:
```bash
# Add custom evaluation to evaluations/all.sh
python /scripts/evaluations/custom_metrics/my_metric.py \
    --input-dir /data/gold \
    --output-file /outputs/metrics/custom_results.csv \
    --methods ranked1,ranked2,combat

# Integrate with existing evaluation framework
bash /scripts/evaluations/classify_batch_bio_within_dataset/classify.sh custom_method
```

## Performance Monitoring and Resource Management

## Troubleshooting and Diagnostics

### Systematic Problem Diagnosis

**1. Pipeline Initialization Issues**:
```bash
# Verify environment setup
echo "=== Environment Check ==="
echo "PYTHONPATH: $PYTHONPATH"
echo "PATH: $PATH"
echo "Working directory: $(pwd)
### Emergency Recovery Procedures

**Pipeline Reset**:
```bash
# Complete pipeline reset
echo "Performing complete pipeline reset..."
rm -rf /data/.cache/*
rm -f /outputs/*.log
find /data/gold -name "*.csv" ! -name "unadjusted.csv" -delete
echo "Reset complete. Restart pipeline with: ./scripts/all.sh"
```

**Partial Recovery**:
```bash
# Recover from specific phase failure
recover_from_adjust_failure() {
    echo "Recovering from adjustment phase failure..."
    rm -rf /data/.cache/gmm_cache/
    rm -f /outputs/adjust.log
    find /data/gold -name "ranked*.csv" -delete
    echo "Ready to restart adjustment phase"
}

# Recover from evaluation failure
recover_from_eval_failure() {
    echo "Recovering from evaluation phase failure..."
    rm -rf /data/.cache/classify_hashes/
    rm -f /outputs/metrics/*.csv
    rm -f /outputs/metrics.log
    echo "Ready to restart evaluation phase"
}
```

---

> **Navigation**: [← Main README](../README.md) | [Batch Correction Methods →](adjust/README.md) | [Evaluation Framework →](evaluations/README.md) | [Data Preparation →](prepdata/README.md)