# Evaluation Framework

> **Navigation**: [← Main README](../../README.md) | [← Pipeline Documentation](../README.md) | [← Batch Correction Methods](../adjust/README.md) | [Data Preparation →](../prepdata/README.md)

Evaluation framework for assessing batch correction methods through classification, statistical metrics, and specialized analyses.

## Classification-Based Evaluation

**Intra-Dataset** (`classify_batch_bio_within_dataset/`) - Batch vs. biological signal classification
**Cross-Dataset** (`classify_er_mixed_datasets/`) - ER status prediction across studies
**Feature Importance** (`classifier_feature_importance/`) - Feature stability analysis

## Statistical Metrics (`small_evals/`)

**MMD** (`mmd.py`) - Distribution similarity via kernel methods (lower is better: <0.05 excellent, 0.05-0.15 good, >0.25 poor)
**MSE** (`mse.py`) - Reconstruction quality (balance metric for over/under-correction)
**Mutual Information** (`mutual_info.py`) - Batch-biology entanglement (Batch-Biology MI <0.1 excellent, preservation ratio >0.8)

## Specialized Analyses

**BatchQC** (`batchqc/`) - Quality control with PCA, heatmaps, statistical testing
**ESR1 Analysis** (`esr1/`) - Estrogen receptor pathway preservation
**Dimensionality Reduction** (`reduce_data_for_viewing/`) - t-SNE, UMAP, PCA visualization

## Robustness Testing (`robustifying/`)

**Simulation Studies** - Controlled batch effects with parameter sweeps across platforms
**Real Data Validation** - TB and cancer studies with cross-study prediction
**Key Files**: `1_simpipe.R` (simulation), `2_TB_getdata.R` (data prep), `3-6_real_data_pipe.R` (validation)

## Usage

```bash
# Run all evaluations
./scripts/evaluations/all.sh

# Specific modules
./scripts/evaluations/classify_batch_bio_within_dataset/classify.sh
python classify_batch_bio_within_dataset/classify.py --method autoclass --dataset gse19615_gse20194
Rscript robustifying/code/3_real_data_pipe.R

# SLURM submission
./run_in_apptainer.sh --sbatch scripts/evaluations/robustifying/code/3_real_data_pipe.R
```

## Outputs

**Metrics** (`outputs/metrics/`) - CSV files for batch/biological classification, MI, MSE, MMD, ER classification
**Figures** (`outputs/figures/`) - Classification plots, MSE/MMD comparisons, dimensionality reduction
**Tables** (`outputs/tables/`) - LaTeX formatted results


---

> **Navigation**: [← Main README](../../README.md) | [← Pipeline Documentation](../README.md) | [← Batch Correction Methods](../adjust/README.md) | [Data Preparation →](../prepdata/README.md)