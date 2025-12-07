# Pipeline Documentation

> **Navigation**: [← Main README](../README.md) | [Batch Correction Methods →](adjust/README.md) | [Evaluation Framework →](evaluations/README.md) | [Data Preparation →](prepdata/README.md)

Analysis pipeline with four main phases: data preparation, batch correction, evaluation, and visualization.

## Pipeline Phases

**1. Data Preparation** (`prepdata/all.sh`) - Download from OSF/Google Drive/Refinebio, format standardization, gene mapping. See [prepdata/README.md](prepdata/README.md)

**2. Batch Correction** (`adjust/all.sh`) - GMM variants, deep learning (AutoClass, ICVAE, VFAE, Wasserstein), statistical methods (ComBat, MNN, Seurat, LIGER, limma), ranking methods. Parallel/sequential execution based on resource needs. See [adjust/README.md](adjust/README.md)

**3. Evaluation** (`evaluations/all.sh`) - Classification metrics, statistical measures (MMD, MSE, MI), BatchQC, pathway analysis, robustness testing. See [evaluations/README.md](evaluations/README.md)

**4. Visualization** - Performance plots, dimensionality reduction (t-SNE, UMAP, PCA), LaTeX tables, confusion matrices

## Execution

**Master Script**: `scripts/all.sh` with error handling, environment setup, logging, modular phase control

```bash
# Docker
./run_docker.sh

# Apptainer
./run_in_apptainer.sh shell
./run_in_apptainer.sh scripts/all.sh
./run_in_apptainer.sh --sbatch --time 04:00:00 --mem 128G scripts/adjust/adjustR_data.sh
```

## Caching

**Components**: HashCache (MD5 validation), DataFrameCache (in-memory), file-based locking
**Locations**: `data/.cache/` (gmm_cache/, gdown/, R/, classify_hashes/, [method]_cache.json)

```bash
# Clear caches
rm -rf data/.cache/                    # All caches
rm -rf data/.cache/gmm_cache/          # Method-specific
rm -rf data/.cache/gdown/              # Downloads only
```

## Logging

**Phase Logs**: prepdata.log, prepdata2.log, adjust.log, metrics.log, figures.log
**Method Logs**: esr1_analysis.log, classify.log, hist_gradient_er.log, feature_importance.log

## Configuration

**Data Preparation** (`prepdata/all.sh`) - Enable/disable download sources, adjust parallelism
**Batch Correction** (`adjust/adjustR_data.sh`) - Select methods (parallel/sequential/target), configure datasets, set batch/target columns
**Evaluation** (`evaluations/all.sh`) - Enable/disable evaluation components

---

> **Navigation**: [← Main README](../README.md) | [Batch Correction Methods →](adjust/README.md) | [Evaluation Framework →](evaluations/README.md) | [Data Preparation →](prepdata/README.md)
