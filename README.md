# Analysis code for "Adversarial deep neural networks remove nonlinear batch effects from gene-expression data"

To execute the analysis, you must install [Docker Desktop](https://www.docker.com/products/docker-desktop). Then invoke the `run_docker.sh` script at the command line.

## Pipeline Overview

The analysis pipeline consists of four main phases executed by `scripts/all.sh`:

1. **Data Preparation** (`scripts/prepdata/all.sh`) - Downloads and processes datasets
2. **Batch Effect Adjustment** (`scripts/adjust/all.sh`) - Applies various adjustment methods  
3. **Metrics Calculation** (`scripts/metrics/all.sh`) - Computes evaluation metrics
4. **Figure Generation** (`scripts/figures/all.sh`) - Creates visualizations and plots

All phases are currently enabled and will run automatically when executing the pipeline.

## Dataset Download System

The project uses a modular dataset download system that separates downloading from OSF and Google Drive sources. The system automatically organizes files into the proper directory structure for downstream processing.

### Docker Usage

#### Full Pipeline (Recommended)
```bash
# Run the complete analysis pipeline (prepdata, adjust, metrics, figures)
./run_docker.sh
```

#### Running Individual Pipeline Phases

Comment out lines in the various all.sh files as needed.
The scripts are designed to be idempotent; running twice will give the same result as running once.
Many scripts are designed cache their results; if the input has not changed, the execution will be skipped, though this is not true for all scripts.
Commenting out lines is therefore required to hasten iteration.

### Pipeline Configuration

The analysis pipeline can be customized by modifying the scripts in each phase:

- **Data Preparation**: Downloads datasets from OSF and Google Drive, organizes files, converts formats, and generates dataset combinations
- **Batch Adjustment**: Applies multiple adjustment methods including AutoClass, ICVAE, VFAE, Wasserstein, and various R-based methods (gmm_adjust uses the R pipeline)
- **Metrics**: Computes mutual information, classification accuracy, MSE, MMD, and other evaluation metrics
- **Figures**: Generates ER classification plots, t-SNE visualizations, and MSE/MMD comparison figures

### Dataset Sources

The system supports downloading from two main sources:

- **OSF (Open Science Framework)**: Primary datasets including GSE19615, GSE20194, GSE20271, GSE23720, GSE25055, GSE25065, GSE31448, GSE45255, GSE58644, GSE62944_Tumor, GSE76275, GSE81538, GSE96058_HiSeq, GSE96058_NextSeq, METABRIC
- **Google Drive**: Additional datasets including GSE115577, GSE123845, GSE163882

#### Debugging

1. **Check pipeline logs**:
   ```bash
   # View individual phase logs
   cat outputs/prepdata.log
   cat outputs/adjust.log
   cat outputs/metrics.log  
   cat outputs/figures.log
   ```

## Pipeline Execution Order

The `scripts/all.sh` script executes the following phases in order:

1. **Data Preparation**
   - Downloads datasets from OSF (GSE19615, GSE20194, GSE20271, GSE23720, GSE25055, GSE25065, GSE31448, GSE45255, GSE58644, GSE62944_Tumor, GSE76275, GSE81538, GSE96058_HiSeq, GSE96058_NextSeq, METABRIC)
   - Downloads datasets from Google Drive (GSE115577, GSE123845, GSE163882)
   - Organizes and converts raw files to standardized formats
   - Generates all possible dataset combinations for analysis

2. **Batch Effect Adjustment**
   - AutoClass: Automated batch effect correction
   - ICVAE: Information-constrained variational autoencoder
   - VFAE: Variational fair autoencoder  
   - Wasserstein: Wasserstein distance-based correction
   - R-based methods: Traditional statistical adjustment approaches

3. **Metrics Calculation**
   - Mutual information analysis between batches and biological variables
   - Classification accuracy for biological outcomes
   - Mean squared error (MSE) and maximum mean discrepancy (MMD)
   - Histogram gradient boosting for ER status prediction

4. **Figure Generation**
   - ER classification performance plots
   - t-SNE dimensionality reduction visualizations
   - MSE/MMD comparison figures
   - Classification accuracy comparisons across methods

## Data Structure

The pipeline organizes data into several key directories:

### `/data/` - Primary Data Storage

```
data/
├── .cache/              # Cached computation results and hashes
│   ├── gdown/          # Google Drive download cache
│   ├── gmm_cache/      # GMM adjustment method cache
│   ├── matplotlib/     # Matplotlib figure cache
│   ├── R/              # R computation cache
│   └── *.hashes.json   # Hash files for dataset combinations and classifications
├── annotations/         # Gene annotation and mapping files
│   ├── entrez_to_symbol_map.csv    # Gene ID mappings
│   └── GPL96-annotation*.csv       # Platform annotation files
├── combined_data/       # Merged datasets for batch effect analysis
│   └── gse*_gse*/      # Pairwise dataset combinations (e.g., gse19697_gse20194/)
├── gold/               # Processed, analysis-ready datasets
│   ├── gse*/           # Individual GSE datasets (e.g., gse20194/, gse24080/)
│   ├── metabric/       # METABRIC breast cancer dataset
│   └── refine/         # Refinebio processed data
├── raw_data/           # Intermediate processed data from raw downloads
│   └── gse*/           # Dataset-specific processed files
├── raw_download/       # Original downloaded files from OSF/Google Drive
│   └── gse*/           # Raw dataset files as downloaded
├── synthetic/          # Synthetic datasets for method validation
│   ├── *_dims_*_bio_*_batch/  # Synthetic data with varying dimensions and batch effects
│   ├── quad2d/         # 2D quadratic synthetic data
│   ├── reduced_data/   # Dimensionally reduced synthetic datasets
│   └── structured_synthetic/   # Structured synthetic data for testing
└── refinebio.h5        # Refinebio HDF5 data file
```

### `/outputs/` - Analysis Results

```
outputs/
├── *.log               # Pipeline execution logs (prepdata, adjust, metrics, figures)
├── performance_log.csv # Performance tracking across runs
├── favorites/          # Curated best results and figures
├── figures/            # Generated visualizations
│   ├── classification/ # ER status classification performance plots
│   ├── classification_global/ # Global classification results
│   ├── combined_global/ # Combined dataset global analysis
│   ├── histograms/     # Data distribution histograms
│   ├── mse_mmd/        # MSE and MMD comparison plots
│   ├── pca/            # Principal component analysis plots
│   ├── prop/           # Proportion and distribution plots
│   ├── reduced/        # Dimensionally reduced data visualizations
│   ├── scatter_plots/  # Scatter plot analyses
│   ├── violin_plots/   # Distribution violin plots
│   └── er_status_model_comparison_*.pdf # Method comparison plots
├── metrics/            # Quantitative evaluation results
│   ├── batch_classification.csv     # Batch effect classification results
│   ├── true_classification.csv      # Biological classification results
│   ├── mutual_info.csv             # Mutual information metrics
│   ├── mse.csv                     # Mean squared error results
│   ├── mmd.csv                     # Maximum mean discrepancy results
│   ├── hist_gradient_er_*.csv      # ER status prediction results
│   └── feature_importance/         # Feature importance analysis
├── optimizations/      # Model hyperparameter optimization results
└── tables/            # LaTeX formatted result tables
    ├── mmd.tex        # MMD results table
    └── mse.tex        # MSE results table
```

## Docker Container Details

The analysis uses a multi-stage Docker build with:

- **Base**: Bioconductor/bioconductor:RELEASE_3_21
- **Python Environment**: Miniforge3 with scientific computing packages
  - NumPy, scikit-learn, pandas, TensorFlow, PyTorch
  - AIF360, tabulate, UMAP, t-SNE, seaborn
  - OSF client and Google Drive tools
- **R Environment**: Comprehensive R package installation
  - Bioconductor packages for genomics analysis
  - Statistical and machine learning packages
- **External Tools**: AutoClass repository

## Batch Effect Adjustment Methods (`/scripts/adjust/`)

The batch effect adjustment phase implements multiple methods for removing unwanted technical variation while preserving biological signal. The directory contains both deep learning and traditional statistical approaches:

### Core Adjustment Scripts

- **`all.sh`** - Master script that executes all adjustment methods in sequence
- **`adjustR_data.sh`** - Runs R-based methods on individual datasets in parallel/sequential modes
- **`adjustR_individual_prep.sh`** - Applies R methods to individual datasets with global batch correction
- **`adjustR_combined_data.sh`** - Processes combined dataset pairs for cross-study batch correction

### Deep Learning Methods

#### AutoClass
- **`autoclass.py`** - Main AutoClass implementation for automated batch effect correction
- **`autoclass.sh`** - Shell wrapper for running AutoClass on specific datasets
- **`invert_autoclass.py`** - Utility functions for AutoClass batch correction and imputation

#### Information-Constrained Variational Autoencoder (ICVAE)
- **`icvae.py`** - Core ICVAE model with Treeish classifier for batch-aware latent representations
- **`icvae.sh`** - Execution script for ICVAE method on multiple datasets
- **`run_icvae.py`** - Command-line interface for ICVAE training and inference

#### Variational Fair Autoencoder (VFAE)
- **`vfae.py`** - VFAE implementation with encoder/decoder architecture for fair representations
- **`vfae.sh`** - Shell script for running VFAE batch correction
- **`run_vfae.py`** - Training and inference pipeline for VFAE

#### Wasserstein Distance-Based Correction
- **`wasserstein.py`** - Wasserstein distance-based batch effect removal using adversarial training
- **`wasserstein.sh`** - Execution wrapper for Wasserstein method

### Traditional Statistical Methods (R-based)

#### Core R Framework
- **`adjust.R`** - Main R adjustment framework supporting multiple statistical methods:
  - **Combat** - Empirical Bayes batch correction
  - **Quantile normalization** - Distribution-based normalization
  - **FastMNN** - Mutual nearest neighbors batch correction
  - **Seurat integration** - Single-cell RNA-seq integration methods
  - **LIGER** - Integrative non-negative matrix factorization
  - **Limma** - Linear modeling for batch correction
  - **FairAdapt** - Fairness-aware adjustment

#### Gaussian Mixture Model (GMM) Methods
- **`gmm_adjust.R`** - High-level interface for GMM-based batch correction with caching
- **`gmm_parameters.R`** - Parameter estimation and model fitting for GMM adjustment
- **`gmm_transforms.R`** - Data transformation utilities for GMM methods

#### Utility Scripts
- **`subset_data.R`** - Data subsetting utilities for testing and validation

### Method Categories

**Parallel Processing Methods** (run simultaneously):
- Quantile normalization (`npn`)
- Min-mean normalization

**Sequential Processing Methods** (run one at a time):
- FastMNN - Memory-intensive mutual nearest neighbors
- LIGER - Integrative matrix factorization
- GMM variants - Gaussian mixture model approaches

**Target-Specific Methods**:
- Combat - Empirical Bayes correction
- FairAdapt - Causal fairness adjustment
- Limma - Linear modeling approaches

Each method produces adjusted datasets that remove batch effects while attempting to preserve biological variation, enabling comparison of adjustment effectiveness.

## Support

Please log an [issue](https://github.com/srp33/confounded_analysis/issues) if you run into a problem with the analysis. You can also contact us [directly](https://biology.byu.edu/piccolo-lab/contact).
