# Important Feature Selection Documentation

## Workflow
After running the prepdata files and creating gold/, this pipeline performs the following steps:
1. Combines all data, keeping all metadata columns
2. Subsets the combined dataset into target datasets, metabric and gse62944
3. Aligns metadata for target datasetes (metabric and gse62944), creating combined columns for er_status, etc.
4. Adjusts the datasets using each of the provided batch effect adjustment methods
5. Trains a classifier model on each adjusted dataset and saves the model and performance
6. Calculates feature importance values, either t-test or permutation importance
7. Plots heatmaps
8. Ranks and selects top genes

This workflow is also visible in the directory's Snakefile.

## Input Files
The gold/ directory and the names of the target datasets are required, as well as which columns will be aligned and what importance test will be used. align_metadata.py must be customized to combine the columns and map values to binary. 

**Bolded Heading**
