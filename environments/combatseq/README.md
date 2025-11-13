# ComBat-seq Environment

R environment for batch correction analysis using ComBat-seq and related methods.

## Configuration

- **R Version**: 4.0.5
- **Bioconductor**: 3.11
- **CRAN Snapshot**: 2020-06-01 (for reproducibility)

## Key Packages

- **Batch Correction**: sva (ComBat), batchelor, RUVSeq, BatchQC
- **Differential Expression**: DESeq2, limma
- **Single-cell**: Seurat
- **Visualization**: ComplexHeatmap, ggplot2, gridExtra

## Usage

Activate the environment:
```bash
source environments/load_envs.sh combatseq
```

## Notes

This environment uses older package versions (circa 2020) for compatibility with legacy analysis code. For new projects, consider using the `book_chapter` environment with more recent packages.

## Packages Not Included

The following packages from the original Apptainer definition are not included due to installation issues or unavailability:
- `fairadapt` - GitHub-only package, may require manual installation
- `rliger` - May have compatibility issues with R 4.0

Install manually if needed:
```r
remotes::install_github("YosefLab/fairadapt")
remotes::install_github("MacoskoLab/liger")
```
