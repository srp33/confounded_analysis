# ER Classification Scaling Experiment

This pipeline evaluates batch effect correction methods across multiple breast cancer datasets using a scaling experiment design.

## Configuration

All paths can be configured via `config.yaml` or command-line overrides:

```yaml
output_folder: "outputs"           # Where results are saved
data_folder: "data"                # Where input data is located
scripts_dir: "."                   # Directory containing pipeline scripts
adjust_script: "../../adjust/adjust.R"  # Path to adjustment functions
metadata_file: "geo_metadata.csv"  # Dataset metadata file
```

## Running the Pipeline

### Default paths (from config.yaml)
```bash
./run_snakemake.sh
```

### Custom paths via command line
```bash
./run_snakemake.sh --config \
    output_folder=/grphome/grp_batch_effects/outputs \
    data_folder=/grphome/grp_batch_effects/data \
    metadata_file=/path/to/geo_metadata.csv
```

### Running locally (without SLURM)
```bash
snakemake --cores 4 --config output_folder=outputs data_folder=data
```

## Pipeline Steps

1. **combine**: Merges all unadjusted.csv files from gold directory
2. **make_order_files**: Creates randomized training order for each test dataset
3. **subset**: Creates subsets with k training studies + 1 test study
4. **adjust**: Applies batch correction methods to each subset
5. **classify**: Trains/tests ER status classifier on adjusted data
6. **aggregate_metrics**: Combines all classification metrics
7. **plot_performance**: Generates scaling performance plots

## Output Structure

```
{output_folder}/
├── metrics/
│   ├── gmm/
│   ├── min_mean/
│   └── ...
├── classify_er_all/
│   └── all_metrics.csv
└── classify_all_figures/
    ├── absolute_scaling_MCC.png
    ├── absolute_scaling_ROC_AUC.png
    ├── scaling_MCC.png
    └── scaling_ROC_AUC.png
```
