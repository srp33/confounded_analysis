# BatchQC Analysis

This folder contains scripts for comprehensive batch effect analysis and visualization using the BatchQC R package.

## Files

- `batchqc_analysis.R` - Main analysis script that takes input file and output directory as parameters
- `run_batchqc.sh` - Wrapper script that processes all paired datasets
- `test_batchqc.R` - Simple test script for single dataset analysis
- `final_batchqc.R` - Standalone script with hardcoded paths (for testing)

## Usage

### Run analysis on all datasets:
```bash
./run_in_apptainer.sh /scripts/evaluations/batchqc/run_batchqc.sh
```

### Run analysis on a single dataset:
```bash
./run_in_apptainer.sh /scripts/evaluations/batchqc/batchqc_analysis.R /data/paired_datasets/dataset_name/unadjusted.csv /grp_batch_effects/outputs/batchqc/dataset_name_unadjusted
```

### Test with sample data:
```bash
./run_in_apptainer.sh /scripts/evaluations/batchqc/test_batchqc.R
```

## Output

Results are saved to `/grp_batch_effects/outputs/batchqc/` with the following structure:
```
/grp_batch_effects/outputs/batchqc/
├── dataset1_unadjusted/
│   ├── pca_plot.png
│   ├── explained_variation.png
│   ├── sample_correlation_heatmap.png
│   ├── pc1_boxplot.png
│   ├── pc2_boxplot.png
│   └── summary.txt
├── dataset1_gmm_affine/
│   └── ...
└── dataset2_unadjusted/
    └── ...
```

## Generated Visualizations

1. **PCA Plot** - Shows batch separation in principal component space
2. **Explained Variation** - Quantifies variance explained by batch vs biological factors
3. **Sample Correlation Heatmap** - Visualizes sample-to-sample correlations with batch annotations
4. **PC Boxplots** - Shows distribution of principal components by batch

## Requirements

- BatchQC R package (installed in container)
- pheatmap R package (for correlation heatmaps)
- Input data with `meta_source` column for batch information
- Optional `meta_er_status` column for biological condition