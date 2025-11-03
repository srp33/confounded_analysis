# Evaluation Framework Documentation

> **Navigation**: [← Main README](../../README.md) | [← Pipeline Documentation](../README.md) | [← Batch Correction Methods](../adjust/README.md) | [Data Preparation →](../prepdata/README.md)

This document provides documentation for the evaluation framework used to assess batch correction methods. 

## Table of Contents

1. [Evaluation Strategy](#evaluation-strategy)
2. [Classification-Based Evaluation](#classification-based-evaluation)
3. [Statistical Validation Metrics](#statistical-validation-metrics)
4. [Specialized Analysis Modules](#specialized-analysis-modules)
5. [Robustness and Validation Framework](#robustness-and-validation-framework)
7. [Usage Examples](#usage-examples)
8. [Output Structure](#output-structure)

## Evaluation Strategy

### Multi-Dimensional Assessment Philosophy

The evaluation framework employs multiple approaches to assess batch correction effectiveness:

1. **Classification-Based Evaluation**: Quantify batch effect removal and biological signal preservation through supervised learning
2. **Statistical Validation Metrics**: Measure distribution similarity and information preservation using unsupervised approaches
3. **Specialized Analysis Modules**: Domain-specific assessments including pathway analysis and quality control
4. **Robustness Testing**: Method stability, parameter sensitivity, and cross-dataset generalization analysis

### Core Evaluation Principles

- **Dual Objective Assessment**: Balance between batch effect removal and biological signal preservation
- **Method-Agnostic Design**: Evaluation metrics independent of correction method implementation
- **Scalable Framework**: Efficient computation for large-scale comparative studies
- **Reproducible Results**: Deterministic evaluation with proper random seed management
- **Interpretable Outputs**: Clear metrics with biological and statistical significance

## Classification-Based Evaluation

### Intra-Dataset Analysis
- **Directory**: `classify_batch_bio_within_dataset/`
- **Purpose**: Batch vs. biological signal classification within individual studies
- **Files**: `classify.py`, `classify.sh`, `classification_figures.R`
- **Metrics**: Classification accuracy, feature importance, confusion matrices

### Cross-Dataset Prediction
- **Directory**: `classify_er_mixed_datasets/`
- **Purpose**: ER status prediction across different studies using histogram gradient boosting
- **Files**: `HGB_classify_er.py`, `hist_gradient_er_classification.py`, `confusion_matrix.py`
- **Metrics**: Cross-study prediction accuracy, model transferability

### Feature Importance Analysis
- **Directory**: `classifier_feature_importance/`
- **Purpose**: Stability analysis of predictive features across correction methods
- **Files**: `feature_importance_analysis.py`, `feature_importance.sh`
- **Output**: Feature stability rankings, method comparison matrices

## Statistical Validation Metrics

**Directory**: `small_evals/`
**Core Philosophy**: Unsupervised assessment of batch correction quality through distribution analysis

### Maximum Mean Discrepancy (MMD)
- **File**: `mmd.py`
- **Purpose**: Measure distribution similarity between batches after correction using kernel methods
- **Mathematical Foundation**: Two-sample test statistic based on reproducing kernel Hilbert spaces
- **Kernel Strategy**: Multi-scale Gaussian kernels (σ = m/2, m, 2m where m is median distance)
- **Interpretation**: 
  - Range: [0, ∞), Lower values indicate better batch mixing
  - Excellent: < 0.05, Good: 0.05-0.15, Poor: > 0.25
- **Usage**: 
  ```bash
  python small_evals/mmd.py --input-dir data/paired_datasets/gse19615_gse20194/ \
                            --batch-col meta_batch --output-path outputs/metrics/mmd.csv
  ```
- **Output**: Per-dataset MMD scores with method comparisons

### Mean Squared Error (MSE)
- **File**: `mse.py`
- **Purpose**: Quantify data reconstruction quality and information preservation
- **Calculation**: MSE between original and corrected expression matrices
- **Interpretation**: 
  - Balance metric: Too low may indicate over-correction, too high indicates under-correction
  - Optimal range depends on dataset characteristics and noise levels
- **Usage**:
  ```bash
  python small_evals/mse.py --method_comparison --dataset_filter "gse.*" \
                            --output outputs/metrics/mse_comparison.csv
  ```
- **Applications**: Method parameter tuning, correction strength assessment

### Mutual Information Analysis
- **File**: `mutual_info.py`
- **Purpose**: Information-theoretic quantification of batch-biology entanglement
- **Variants Implemented**:
  - **Shannon MI**: Standard mutual information in bits
  - **Probability-based MI**: Continuous prediction compatibility
  - **Determinant-based MI**: Joint probability distribution analysis
- **Key Metrics**:
  - **Batch-Biology MI**: Lower indicates reduced confounding (target: minimize)
  - **Biology-Corrected MI**: Higher indicates preserved signal (target: maximize)
- **Interpretation Guidelines**:
  - Batch-Biology MI: < 0.1 bits (excellent), 0.1-0.3 (good), > 0.5 (poor)
  - Signal preservation ratio: (Biology-Corrected MI) / (Original Biology MI) > 0.8
- **Usage**:
  ```bash
  python small_evals/mutual_info.py --compute_all --save_results \
                                    --output_dir outputs/metrics/mutual_info/
  ```

## Specialized Analysis Modules

### BatchQC Integration
- **Directory**: `batchqc/`
- **Files**: `batchqc_analysis.R`, `final_batchqc.R`, `run_batchqc.sh`
- **Purpose**: Comprehensive quality control using the BatchQC R package
- **Features**:
  - Principal component analysis of batch effects
  - Heatmap visualization of batch-sample relationships
  - Statistical testing for batch effect significance

### Pathway-Specific Analysis

**ESR1 Analysis**:
- **Directory**: `esr1/`
- **File**: `esr1_analysis.py`
- **Purpose**: Estrogen receptor pathway preservation assessment
- **Features**: Pathway enrichment analysis, gene set preservation

### Dimensionality Reduction and Visualization
- **Directory**: `reduce_data_for_viewing/`
- **Files**: `reduce.py`, `plot_reduced.R`
- **Methods**: t-SNE, UMAP, PCA
- **Purpose**: Visual assessment of batch effect removal and biological structure preservation

## Robustness and Validation Framework

### Comprehensive Simulation Studies
- **Directory**: `robustifying/`
- **Purpose**: Systematic evaluation of method robustness across diverse scenarios
- **Design Philosophy**: Controlled batch effect introduction with known ground truth

**Simulation Framework Components**:
- **Batch Effect Modeling**: Realistic batch effect simulation based on real data characteristics
- **Parameter Sweeps**: Systematic variation of batch sizes, effect magnitudes, and data dimensions
- **Cross-Modal Validation**: Testing across RNA-seq, microarray, and single-cell platforms
- **Noise Robustness**: Performance assessment under varying noise levels

**Key Simulation Files**:
- `code/1_simpipe.R`: Core simulation pipeline with configurable parameters
- `code/qsub_simpipe.qsub`: HPC job submission for large-scale simulation studies
- `code/helper.R`: Utility functions for simulation setup and analysis

### Real-World Validation Studies
- **Directory**: `robustifying/`
- **Datasets**: Tuberculosis (TB) and cancer studies with known biological structure
- **Validation Strategy**: Cross-study prediction and biological pathway preservation

**Real Data Pipeline**:
- `code/2_TB_getdata.R`: Tuberculosis dataset acquisition and preprocessing
- `code/3_real_data_pipe.R`: Comprehensive real data validation workflow
- `code/4_real_data_pipe.R`: Extended validation with additional datasets
- `code/5_real_data_pipe.R`: Cross-platform validation studies
- `code/6_real_data_pipe.R`: Final validation and method comparison

**Biological Validation**:
- **Pathway Preservation**: Gene set enrichment analysis post-correction
- **Differential Expression Consistency**: DE gene overlap across batches
- **Cell Type Identification**: Preservation of cellular identity markers

## Usage Examples

### Running Complete Evaluation
```bash
# Run all evaluation metrics
./scripts/evaluations/all.sh

# Run specific evaluation modules
./scripts/evaluations/classify_batch_bio_within_dataset/classify.sh
./scripts/evaluations/small_evals/mmd.sh
./scripts/evaluations/batchqc/run_batchqc.sh
```

### Custom Evaluation Workflows
```bash
# Evaluate specific method
python classify_batch_bio_within_dataset/classify.py --method autoclass --dataset gse19615_gse20194

# Generate comparison plots
Rscript classify_er_mixed_datasets/combo_er_classification_plots.R

# Run robustness testing
Rscript robustifying/code/3_real_data_pipe.R
```

### Batch Processing
```bash
# Submit evaluation jobs to SLURM
./run_in_apptainer.sh --sbatch scripts/evaluations/robustifying/code/3_real_data_pipe.R

```

## Output Structure

### Metrics Directory (`outputs/metrics/`)
- `batch_classification.csv`: Batch effect classification results
- `true_classification.csv`: Biological classification results
- `mutual_info.csv`: Mutual information metrics
- `mse.csv`: Mean squared error results
- `mmd.csv`: Maximum mean discrepancy results
- `er_classification_*.csv`: ER classification performance metrics

### Figures Directory (`outputs/figures/`)
- `classification/`: ER status classification performance plots
- `mse_mmd/`: MSE and MMD comparison plots
- `reduced/`: Dimensionally reduced data visualizations
- `er_status_model_comparison_*.pdf`: Method comparison plots

### Tables Directory (`outputs/tables/`)
- `mmd.tex`: MMD results table (LaTeX format)
- `mse.tex`: MSE results table (LaTeX format)


---

> **Navigation**: [← Main README](../../README.md) | [← Pipeline Documentation](../README.md) | [← Batch Correction Methods](../adjust/README.md) | [Data Preparation →](../prepdata/README.md)