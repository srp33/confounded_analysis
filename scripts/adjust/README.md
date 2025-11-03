# Batch Correction Methods

> **Navigation**: [← Main README](../../README.md) | [← Pipeline Documentation](../README.md) | [Evaluation Framework →](../evaluations/README.md) | [Data Preparation →](../prepdata/README.md)

Documentation for batch effect correction methods available in the pipeline, with particular focus on the gmm_adjust method and its variants, alongside deep learning approaches and statistical methods.

## Method Categories

### Gaussian Mixture Model Suite (Primary Methods)

#### Standard GMM (gmm_adjust)
- **Files**: `gmm_adjust.R`, `gmm_adjust_python.py`
- **Description**: 1D 2-component Gaussian Mixture Model with posterior-mean priors - the primary method for batch effect correction in this pipeline
- **Architecture**: Expectation-Maximization algorithm with Bayesian priors
- **Key Features**: 
  - Dirichlet priors for mixture weights (configurable alpha >= 1.0)
  - Inverse-Gamma priors for variances (alpha > 1 for posterior mean existence)
  - Cross-platform compatibility (Python and R implementations)
  - Result caching with hash-based validation
  - Rigorous testing across multiple datasets and evaluation metrics
- **Parameters**:
  - `max_iter`: Maximum EM iterations (default: 100)
  - `tol`: Convergence tolerance for log-likelihood (default: 1e-4)
  - `weight_alpha`: Dirichlet prior pseudo-count for mixture weights
  - `variance_alpha`: Inverse-Gamma prior shape parameter for variances

#### Nonlinear GMM
- **File**: `gmm_adjust_nonlinear.R`
- **Description**: Extended variants for complex, non-linear batch structures
- **Features**: 
  - Handles complex batch effect patterns
  - Non-linear transformation capabilities
  - Mixture modeling for heterogeneous data

#### Global GMM
- **File**: `gmm_global_simple.R`
- **Description**: Simplified approach for large-scale corrections
- **Features**: 
  - Optimized for computational efficiency
  - Reduced memory footprint
  - Supports large datasets
  - Streamlined parameter estimation

### Deep Learning Approaches

#### AutoClass Framework
- **Files**: `autoclass.py`, `autoclass.sh`, `invert_autoclass.py`
- **Description**: Neural networks with self-configuring architecture for batch effect removal using adversarial training
- **Architecture**: Encoder-decoder with adversarial discriminator for batch effect removal
- **Features**: 
  - Data normalization detection (microarray vs RNA-seq)
  - Adversarial training with configurable weight (default: 0.002)
  - Encoder layer customization (default: [128])
  - Data imputation with cellwise normalization options
- **Parameters**:
  - `encoder_layer_size`: Neural network architecture (default: [128])
  - `adversarial_weight`: Adversarial loss weight (default: 0.002)
  - `epochs`: Training epochs (default: 400)
  - `lr`: Learning rate (default: 15)
  - `reg`: Regularization strength (default: 0.0001)
  - `dropout_rate`: Dropout probability (default: 0.2)

#### Variational Autoencoders

**ICVAE (Information-Constrained VAE)**
- **Files**: `icvae.py`, `icvae.sh`, `run_icvae.py`
- **Description**: Information-constrained VAE with Treeish classifier for batch-aware latent representations
- **Architecture**: VAE with auxiliary classifier and mutual information penalty
- **Features**: 
  - Treeish classifier mimicking decision tree behavior for batch prediction
  - Controlled information flow through mutual information penalty
  - Feature masking for random feature selection
  - Biological signal preservation through information constraints
- **Parameters**:
  - `--latent-dim`: Latent space dimensionality (default: 10)
  - `--hidden-dim`: Hidden layer size for VAE (default: 128)
  - `--hidden-dim-aux`: Hidden layer size for auxiliary classifier (default: 64)
  - `--epochs`: Training epochs (default: 100)
  - `--learning-rate`: Optimizer learning rate (default: 1e-3)
  - `--batch-size`: Training batch size (default: 64)
  - `--w-kl`: KL divergence loss weight (default: 1.0)
  - `--w-mi-penalty`: Mutual information penalty weight (default: 1.0)

**VFAE (Variational Fair Autoencoder)**
- **Files**: `vfae.py`, `vfae.sh`, `run_vfae.py`
- **Description**: Ensures equitable correction across biological subgroups using Maximum Mean Discrepancy (MMD) penalty
- **Architecture**: VAE with fairness constraints through MMD loss
- **Features**: 
  - Fairness constraints using MMD penalty
  - RBF kernel-based distribution matching
  - Bias mitigation across sensitive attributes
  - Equitable representation learning
- **Parameters**:
  - `--latent-dim`: Latent space dimensionality (default: 10)
  - `--hidden-dim`: Hidden layer size (default: 128)
  - `--epochs`: Training epochs (default: 100)
  - `--learning-rate`: Optimizer learning rate (default: 1e-3)
  - `--batch-size`: Training batch size (default: 64)
  - `--w-kl`: KL divergence loss weight (default: 1.0)
  - `--w-mmd`: MMD penalty weight (default: 10.0)
  - `--mmd-gamma`: RBF kernel gamma parameter (default: 1.0)

**Wasserstein Adversarial Training**
- **Files**: `wasserstein.py`, `wasserstein.sh`
- **Description**: Adversarial training using Wasserstein distance for stable batch correction
- **Architecture**: Generator-critic architecture with Wasserstein loss
- **Features**: 
  - Wasserstein distance-based adversarial training
  - Gradient penalty for training stability
  - Multi-layer perceptron critic network
  - Less sensitive to hyperparameter choices

### Statistical and Traditional Methods

#### R-Based Methods Integration
- **Files**: `adjust.R`, `adjustR_data.sh`, `adjustR_individual_prep.sh`, `adjustR_paired_datasets.sh`
- **Description**: Comprehensive R-based statistical methods with parallel and sequential execution modes
- **Execution Modes**:
  - **Parallel Methods**: `ranked1`, `ranked2`, `ranked_batch` (all datasets processed simultaneously)
  - **Sequential Methods**: `mnn`, `liger` (one dataset at a time)
  - **Target Methods**: `combat`, `fairadapt`, `limma` (specialized execution)

#### Implemented Statistical Methods

**ComBat (Empirical Bayes)**
- **Implementation**: Via `sva` package in R
- **Description**: Gold standard empirical Bayes batch correction method
- **Features**: Works with small sample sizes, preserves biological variation
- **Use Case**: General-purpose batch correction for microarray and RNA-seq data

**Quantile Normalization**
- **Implementation**: Via `preprocessCore` package
- **Description**: Distribution-based normalization for cross-platform studies
- **Features**: Forces identical distributions across batches
- **Use Case**: Cross-platform data integration

**MNN (Mutual Nearest Neighbors)**
- **Implementation**: Via `batchelor` package
- **Description**: Identifies mutual nearest neighbors for batch integration
- **Features**: Preserves local neighborhood structure, handles non-linear batch effects
- **Use Case**: Single-cell data integration, complex batch structures

**Seurat Integration**
- **Implementation**: Via `Seurat` package
- **Description**: Multi-modal integration methods
- **Features**: Canonical correlation analysis, anchor-based integration
- **Use Case**: Single-cell multi-modal data integration

**LIGER (Integrative NMF)**
- **Implementation**: Via `rliger` package
- **Description**: Integrative non-negative matrix factorization
- **Features**: Shared and dataset-specific factors, supports large datasets
- **Use Case**: Large-scale multi-dataset integration

**Limma (Linear Modeling)**
- **Implementation**: Via `limma` package
- **Description**: Linear modeling with explicit batch covariates
- **Features**: Covariate modeling, differential expression integration
- **Use Case**: Controlled batch correction with known covariates

**FairAdapt**
- **Implementation**: Via `fairadapt` package
- **Description**: Causal fairness-aware batch correction
- **Features**: Causal graph-based approach, fairness constraints
- **Use Case**: Bias-aware batch correction with fairness considerations

**Ranked Methods**
- **Implementation**: Custom ranking-based normalization
- **Variants**: 
  - `ranked1`: Basic rank normalization
  - `ranked2`: Enhanced rank normalization
  - `ranked_batch`: Batch-aware rank normalization
- **Features**: Non-parametric, distribution-free approach
- **Use Case**: Normalization for non-normal distributions



## Method Selection and Usage

### Resource Requirements and Execution Strategy

**Memory-Intensive Methods** (Sequential Execution):
- **AutoClass**: GPU-accelerated, high memory usage, 400 epochs default
- **ICVAE**: Moderate memory, auxiliary classifier overhead
- **VFAE**: Moderate memory, MMD computation overhead  
- **Wasserstein**: High memory for critic network, gradient penalty computation

**Lightweight Methods** (Parallel Execution):
- **GMM variants (Recommended)**: Low memory, efficient EM algorithm, rigorously tested
- **Ranked methods**: Minimal memory, distribution-free
- **Statistical methods**: Variable memory based on method complexity

**Execution Control**:
- Parallel methods defined in `ADJUSTERS_PARALLEL` array
- Sequential methods defined in `ADJUSTERS_SEQUENTIAL` array
- Target methods defined in `ADJUSTERS_TARGET` array

### Detailed Usage Examples

#### Deep Learning Methods

```bash
# AutoClass with custom parameters
python /scripts/adjust/autoclass.py \
  -i /data/gold/gse49711/unadjusted.csv \
  -o /data/gold/gse49711/autoclass.csv \
  -b meta_Sex

# ICVAE with custom configuration
python /scripts/adjust/run_icvae.py \
  -i /data/gold/gse24080/unadjusted.csv \
  -o /data/gold/gse24080/icvae.csv \
  -b meta_batch \
  --latent-dim 20 \
  --epochs 400 \
  --w-mi-penalty 1.5

# VFAE with MMD penalty tuning
python /scripts/adjust/run_vfae.py \
  -i /data/gold/gse49711/unadjusted.csv \
  -o /data/gold/gse49711/vfae.csv \
  -b meta_Sex \
  --w-mmd 10.0 \
  --mmd-gamma 1.0

# Wasserstein adversarial training
python /scripts/adjust/wasserstein.py \
  --input /data/gold/dataset/unadjusted.csv \
  --output /data/gold/dataset/wasserstein.csv \
  --batch-col meta_batch
```

#### Statistical Methods

```bash
# Run all parallel R methods
bash /scripts/adjust/adjustR_data.sh

# Run specific R method with dataset
Rscript /scripts/adjust/adjust.R \
  --adjuster combat \
  --dataset gse49711 \
  --input-file /data/gold/gse49711/unadjusted.csv \
  --output-file /data/gold/gse49711/combat.csv

# Run GMM adjustment
Rscript /scripts/adjust/gmm_adjust.R \
  --input /data/gold/dataset/unadjusted.csv \
  --output /data/gold/dataset/gmm.csv \
  --max_iter 100 \
  --tol 1e-4
```

### Method Configuration Details

#### AutoClass Configuration
```python
# Core AutoClass parameters
BatchCorrectImpute(
    genes,                    # Gene expression matrix
    batches,                  # Batch labels
    cellwise_norm=False,      # Cell-wise normalization
    log1p=False,             # Log1p transformation
    verbose=True,            # Verbose output
    encoder_layer_size=[128], # Encoder architecture
    adversarial_weight=0.002, # Adversarial loss weight
    epochs=400,              # Training epochs
    lr=15,                   # Learning rate
    reg=0.0001,              # Regularization strength
    dropout_rate=0.2         # Dropout probability
)
```

#### ICVAE Configuration
```python
# ICVAE hyperparameters
--latent-dim 10              # Latent space dimensionality
--hidden-dim 128             # VAE hidden layer size
--hidden-dim-aux 64          # Auxiliary classifier hidden size
--epochs 100                 # Training epochs
--learning-rate 1e-3         # Optimizer learning rate
--batch-size 64              # Training batch size
--w-kl 1.0                   # KL divergence weight
--w-mi-penalty 1.0           # Mutual information penalty weight
```

#### VFAE Configuration
```python
# VFAE fairness parameters
--latent-dim 10              # Latent space dimensionality
--hidden-dim 128             # Hidden layer size
--epochs 100                 # Training epochs
--learning-rate 1e-3         # Optimizer learning rate
--batch-size 64              # Training batch size
--w-kl 1.0                   # KL divergence weight
--w-mmd 10.0                 # MMD penalty weight
--mmd-gamma 1.0              # RBF kernel gamma parameter
```

#### GMM Configuration
```r
# GMM Bayesian parameters
GaussianMixture1D(
  max_iter = 100,            # Maximum EM iterations
  tol = 1e-4,                # Convergence tolerance
  weight_alpha = NULL,       # Dirichlet prior for weights
  variance_alpha = NULL      # Inverse-Gamma prior for variances
)
```

#### R Statistical Methods Configuration
```r
# Environment configuration
Sys.setenv(OMP_NUM_THREADS = 1)  # Single-threaded OpenMP
ulimit -s unlimited               # Unlimited stack size

# Method-specific parameters
ADJUSTERS_PARALLEL <- c("ranked1", "ranked2", "ranked_batch")
ADJUSTERS_SEQUENTIAL <- c("mnn", "liger")  
ADJUSTERS_TARGET <- c("combat", "fairadapt", "limma")
```


## Performance Optimization

### Scheduling
- **Resource-Optimized Execution**: Parallel processing for lightweight methods, sequential for memory-intensive approaches
- **Adaptive Caching**: Hash-based result storage to avoid redundant computation
- **Performance Monitoring**: Execution time and resource usage tracking

### Caching System
- **Location**: `data/.cache/gmm_cache/`, `data/.cache/R/`
- **Hash-Based Validation**: Content-based caching prevents redundant computations
- **Method-Specific Caches**: Separate caches for different method types

## Method Validation and Quality Assessment

### Validation Framework

#### Synthetic Data Validation
- **Controlled Environments**: Known batch effect characteristics for ground truth comparison
- **Simulation Parameters**: Configurable batch effect strength, biological signal levels
- **Validation Metrics**: Recovery of known biological signals, batch effect removal efficiency
- **Implementation**: Available in `../evaluations/` directory

#### Real Data Cross-Validation
- **Multi-Study Validation**: Cross-study validation using multiple datasets (GSE49711, GSE20194, GSE24080)
- **Hold-Out Testing**: Reserve datasets for independent validation
- **Cross-Platform Assessment**: Microarray vs RNA-seq performance comparison
- **Biological Consistency**: Pathway enrichment and differential expression validation

#### Performance Benchmarking
- **Computational Metrics**: Runtime, memory usage, scalability assessment
- **Resource Profiling**: CPU/GPU utilization, memory peak usage
- **Scalability Testing**: Performance across different dataset sizes
- **Efficiency Comparison**: Method-to-method computational cost analysis

### Quality Assessment Metrics

#### Batch Effect Removal Assessment
```bash
# Classification-based metrics
python /scripts/evaluations/classify_batch_bio_within_dataset/classify.py

# Mutual information analysis
python /scripts/evaluations/small_evals/mutual_info.py

# MMD (Maximum Mean Discrepancy) evaluation
python /scripts/evaluations/small_evals/mmd.py
```

**Key Metrics**:
- **Batch Classification Accuracy**: Lower is better (indicates successful batch removal)
- **Mutual Information**: Measures dependence between corrected data and batch labels
- **MMD Score**: Distributional distance between batches after correction
- **Silhouette Score**: Batch separation in corrected data

#### Biological Signal Preservation Assessment
```bash
# Biological classification accuracy
python /scripts/evaluations/classify_er_mixed_datasets/classify.sh

# Feature importance analysis
python /scripts/evaluations/classifier_feature_importance/feature_importance_analysis.py

# MSE analysis for signal preservation
python /scripts/evaluations/small_evals/mse.py
```

**Key Metrics**:
- **Biological Classification Accuracy**: Higher is better (preserved biological signal)
- **Pathway Enrichment Consistency**: Maintained biological pathways after correction
- **Differential Expression Concordance**: Agreement with known biological differences
- **Feature Importance Preservation**: Retention of biologically relevant features


## Adding New Methods

### Implementation Framework
1. **Method Implementation**: Develop in Python (`*.py`) or R (`*.R`)
2. **Execution Wrapper**: Create shell script (`*.sh`) for resource management
3. **Pipeline Integration**: Update `all.sh` with method orchestration
4. **Evaluation Integration**: Add method assessment in `../evaluations/`
5. **Documentation**: Update method descriptions and usage examples

---

> **Navigation**: [← Main README](../../README.md) | [← Pipeline Documentation](../README.md) | [Evaluation Framework →](../evaluations/README.md) | [Data Preparation →](../prepdata/README.md)