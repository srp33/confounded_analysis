# Data Preparation Pipeline Documentation

> **Navigation**: [← Main README](../../README.md) | [← Pipeline Documentation](../README.md) | [← Batch Correction Methods](../adjust/README.md) | [← Evaluation Framework](../evaluations/README.md)

This document provides documentation for the data preparation pipeline, covering dataset acquisition, processing workflows, data organization strategies, and customization options for adding new data sources.

## Data Preparation Overview

The data preparation pipeline handles the workflow from raw data acquisition to analysis-ready datasets, supporting multiple data sources and formats with intelligent caching, quality control, and processing.

### Pipeline Architecture

The data preparation system follows a modular architecture with five core phases:

1. **Dataset Download**: Multi-source acquisition from OSF, Google Drive, and Refinebio with resume capability
2. **Format Conversion**: Standardization to HDF5, CSV, and RDS formats with gene annotation
3. **Quality Control**: Comprehensive validation including sample filtering and batch detection
4. **Dataset Combination**: Automated pairwise combinations for cross-study batch effect analysis
5. **Synthetic Data Generation**: Controlled datasets with known batch effects for method validation

### Execution Control

**Main Pipeline Script**: `all.sh`
- Orchestrates all data preparation phases with intelligent dependency management
- Supports selective phase execution and caching optimization
- Integrates with container environments (Docker/Apptainer)

**Individual Phase Scripts**:
- `download_datasets.py` - Multi-source dataset acquisition
- `convert_raw_files.py` - Format standardization and conversion
- `combine_all.py` - Automated dataset combination generation
- `generate_structured_synthetic.py` - Synthetic data creation

## Data Sources and Acquisition

### Supported Data Sources

The pipeline supports multiple data repositories with specialized downloaders for each platform:

**OSF (Open Science Framework) Repository**:
- **Primary Collection**: 13 cancer gene expression studies
- **Cancer Studies**: GSE19615, GSE20194, GSE20271, GSE23720, GSE25055, GSE25065, GSE31448, GSE45255, GSE58644, GSE62944_Tumor, GSE76275, GSE81538, METABRIC
- **Platform Studies**: GSE96058_HiSeq, GSE96058_NextSeq (cross-platform batch effect analysis)
- **Access Method**: Direct API integration with project-based organization

**Google Drive Repository**:
- **Recent Studies**: GSE115577, GSE123845, GSE163882
- **Large-Scale Datasets**: TCGA subsets and specialized collections
- **Access Method**: Folder-based organization with authentication support

**Refinebio Integration**:
- **Processed Datasets**: Pre-normalized and quality-controlled expression data
- **Standardization**: Consistent gene annotation and sample metadata
- **Integration**: Automated download and format conversion pipeline

### Download System Architecture

**Modular Downloader Framework**:
- **Base Class**: `downloaders/base.py` - Common interface for all downloaders
- **OSF Downloader**: `downloaders/osf_downloader.py` - OSF API integration
- **Google Drive Downloader**: `downloaders/gdrive_downloader.py` - Drive API with authentication
- **Configuration**: `config.py` - Centralized source configuration and parameters

### Download Procedures

**Complete Dataset Acquisition**:
```bash
# Download all configured datasets
python download_datasets.py --source all

# Download from specific source with project ID
python download_datasets.py --source osf --project-id eky3p

# Download from Google Drive with folder ID
python download_datasets.py --source gdrive --folder-id 1smhpktMRyP4yyFHKHSisxRd9jwb8kvrq

# Download specific datasets only
python download_datasets.py --source osf --datasets GSE20194,GSE24080,METABRIC
```

**Advanced Download Options**:
```bash
# Resume interrupted downloads
python download_datasets.py --resume --source all

# Download with custom configuration
python download_datasets.py --source osf --max-retries 5 --timeout 120

# Verify existing downloads
python download_datasets.py --verify --source all

# Download and organize in single step
python download_datasets.py --source osf --organize
```

**Download Status and Monitoring**:
```bash
# Check download status
python download_datasets.py --status

# Monitor download progress
tail -f outputs/prepdata.log

# Verify download integrity
python download_datasets.py --verify-integrity --source all
```

## Data Processing Pipeline

### Format Conversion and Standardization

**Multi-Format Conversion System**:
- **Primary Converter**: `convert_raw_files.py` - Main conversion orchestrator
- **HDF5 Conversion**: `convert_to_h5.py` - Efficient binary storage format
- **R Integration**: `create_smaller_csv.R`, `make_rds_data.sh` - R-compatible formats
- **Matrix Operations**: `transpose_matrix.py` - Data structure optimization


### Data Organization and File Management

**Directory Structure Management**:
- **Organization Script**: `organize_downloaded_files.py` - Automated file organization
- **Structure Validation**: Consistent directory hierarchy across all datasets
- **File Classification**: Automatic detection of expression, metadata, and annotation files
- **Duplicate Handling**: Intelligent duplicate detection and resolution

**Processing Workflow Integration**:
```bash
# Complete processing pipeline
python convert_raw_files.py --input data/raw_download/ --output data/gold/ --full-pipeline

# Selective processing with quality control
python convert_raw_files.py --dataset gse20194 --quality-check --gene-annotation

# Batch processing with parallel execution
python convert_raw_files.py --parallel --max-workers 4 --datasets gse20194,gse24080

# Processing with custom parameters
python convert_raw_files.py --min-samples 20 --min-genes 1000 --expression-threshold 0.1
```

### Data Flow and Directory Structure

**Processing Pipeline Stages**:
1. **Raw Downloads** (`data/raw_download/`) - Original files from data sources
2. **Organized Raw** (`data/raw_data/`) - Structured and validated raw files
3. **Analysis-Ready** (`data/gold/`) - Processed, standardized datasets
4. **Combined Datasets** (`data/paired_datasets/`) - Cross-study combinations
5. **Synthetic Data** (`data/synthetic/`) - Generated validation datasets

**File Organization Workflow**:
```bash
# Automated file organization from downloads
python organize_downloaded_files.py --source-dir data/raw_download/ --target-dir data/raw_data/

# Create analysis-ready datasets with full processing
python convert_raw_files.py --input data/raw_data/ --output data/gold/ --full-pipeline

# Generate structured directory hierarchy
python organize_downloaded_files.py --create-structure --validate-hierarchy
```

## Dataset Combination Strategy

### Pairwise Dataset Generation

**Combination Generation**:
- **Files**: `combine_datasets.py`, `combine_all.py`, `generate_all_combinations.py`
- **Purpose**: Create all possible two-study combinations for cross-batch validation
- **Features**: Metadata harmonization, sample size balancing, platform integration

**Usage Examples**:
```bash
# Generate all pairwise combinations
python generate_all_combinations.py --input data/gold/ --output data/paired_datasets/

# Combine specific datasets
python combine_datasets.py --dataset1 gse20194 --dataset2 gse24080 --output data/paired_datasets/

# Generate combinations with filtering
python combine_all.py --min_samples 50 --max_combinations 100
```

## Configuration and Customization

### Dataset Configuration

**Configuration File**:
- **File**: `config.py`
- **Purpose**: Central configuration for data sources, processing parameters
- **Features**: Download URLs, processing options, quality thresholds

**Customization Options**:
```python
# config.py example
DOWNLOAD_SOURCES = {
    'osf': ['GSE20194', 'GSE24080', 'METABRIC'],
    'gdrive': ['GSE115577', 'GSE123845'],
    'refinebio': True
}

PROCESSING_OPTIONS = {
    'min_samples': 20,
    'min_genes': 1000,
    'quality_threshold': 0.8
}
```

## Performance Optimization

### Caching System

**Intelligent Caching**:
- **Location**: `data/.cache/gdown/`, `data/.cache/R/`
- **Features**: Hash-based validation, resume capability, parallel processing
- **Benefits**: Avoid redundant downloads, faster processing, robust error recovery

**Cache Management**:
```bash
# Check cache status
du -sh data/.cache/*

# Clear download cache
rm -rf data/.cache/gdown/
```

## Data Structure and Organization

### Input Data Structure
```
data/raw_download/          # Original downloaded files
├── gse20194/              # Individual dataset directories
│   ├── expression_data.csv
│   ├── metadata.csv
│   └── platform_info.txt
└── metabric/
    ├── clinical_data.csv
    └── expression_matrix.csv
```

### Processed Data Structure
```
data/gold/                  # Analysis-ready datasets
├── gse20194/
│   ├── expression.csv     # Standardized expression matrix
│   ├── metadata.csv       # Harmonized sample annotations
│   └── gene_info.csv      # Gene annotation information
└── paired_datasets/        # Combined datasets
    └── gse20194_gse24080/
        ├── combined_expression.csv
        ├── combined_metadata.csv
        └── batch_info.csv
```


## Data Management and Customization

### Adding New Data Sources

**Implementation Framework**:

1. **Create Custom Downloader**:
   ```python
   # Create new file: downloaders/custom_downloader.py
   from .base import BaseDownloader
   
   class CustomDownloader(BaseDownloader):
       def list_available_files(self):
           # Implement source-specific file listing
           pass
       
       def find_dataset_files(self, dataset_id):
           # Implement dataset-specific file discovery
           pass
   ```

2. **Update Configuration**:
   ```python
   # Add to config.py
   CUSTOM_SOURCE_CONFIG = {
       'api_endpoint': 'https://api.custom-source.org',
       'authentication': 'token_based',
       'rate_limit': 100,  # requests per minute
       'supported_formats': ['csv', 'tsv', 'h5']
   }
   ```

3. **Integration with Main Pipeline**:
   ```python
   # Update download_datasets.py
   from downloaders.custom_downloader import CustomDownloader
   
   # Add source option to argument parser
   parser.add_argument('--source', choices=['osf', 'gdrive', 'custom'])
   ```

---

> **Navigation**: [← Main README](../../README.md) | [← Pipeline Documentation](../README.md) | [← Batch Correction Methods](../adjust/README.md) | [← Evaluation Framework](../evaluations/README.md)