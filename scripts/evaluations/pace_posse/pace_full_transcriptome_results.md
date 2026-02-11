# PACE Comprehensive Real Transcriptome Test Results

**Generated:** 2026-01-01 13:37:39  
**Test Script:** test_v34_full_transcriptome.py  
**Dataset:** Real TB transcriptome data from multiple studies

## Executive Summary

**Best Performer:** posse_v50_balanced  
**Version:** POSSE v5.0  
**Method:** ComBat-Initialized  
**Overall Score:** 0.598

🧬 **POSSE v5.0 ComBat-Initialized shows best performance!**

## Detailed Results

| Configuration | Version | Method | Signal↑ | Artifact↑ | Stable↑ | DE↑ | HK↑ | LinSep↑ | Semantic↑ | Score↑ |
|---------------|---------|--------|---------|-----------|---------|-----|-----|---------|-----------|--------|
| posse_v50_balanced | POSSE v5.0 | ComBat-Initialized | -0.218 | 0.019 | 0.013 | 0.705 | 0.784 | 0.947 | 0.617 | 0.598 |
| posse_v50_aggressive | POSSE v5.0 | ComBat-Initialized | -0.218 | 0.019 | 0.013 | 0.705 | 0.760 | 0.947 | 0.617 | 0.597 |
| posse_v50_conservative | POSSE v5.0 | ComBat-Initialized | -0.218 | 0.019 | 0.013 | 0.706 | 0.811 | 0.947 | 0.532 | 0.591 |
| posse_v50_sniper | POSSE v5.0 | ComBat-Initialized | -0.218 | 0.019 | 0.013 | 0.710 | 0.804 | 0.951 | 0.510 | 0.581 |
| v30_activity_focused | PACE v3.0 | Activity-Gated | -0.304 | 0.025 | 0.022 | 0.611 | 0.931 | 0.948 | 0.000 | 0.564 |
| v22_aggressive_centered | PACE v2.2 | Centered Cosine | -0.312 | 0.025 | 0.021 | 0.629 | 0.943 | 0.955 | 0.000 | 0.554 |
| v34_housekeeping_focused | PACE v3.4 | Housekeeping Anchors | -0.295 | 0.022 | 0.017 | 0.597 | 0.946 | 0.952 | 0.000 | 0.552 |
| v35_affine_focused | PACE v3.5 | Affine Anchors | -0.295 | 0.026 | 0.012 | 0.612 | 0.927 | 0.961 | 0.000 | 0.546 |
| pace_v31_pure_local | PACE v3.1 | Pure Local | -0.299 | 0.027 | 0.016 | 0.599 | 0.918 | 0.966 | 0.000 | 0.544 |
| v32_iterative_consensus | PACE v3.2 | Iterative Consensus | -0.311 | 0.026 | 0.014 | 0.614 | 0.903 | 0.942 | 0.000 | 0.544 |
| gmm_baseline | GMM | Gaussian Mixture | -0.335 | 0.014 | 0.007 | 0.702 | 0.816 | 0.990 | 0.137 | 0.521 |
| combat_baseline | ComBat | Standard | -0.327 | 0.018 | 0.008 | 0.411 | 1.000 | 0.500 | 0.000 | 0.459 |

## Metric Interpretation

- **Signal↑**: Signal preservation (-1 to 1, higher better)
  - Measures how well TB vs control biological signal is maintained
  - Based on correlation of corrected data with true biological labels
  - +1.0 = perfect preservation, 0.0 = no correlation, -1.0 = signal inverted

- **Artifact↑**: Artifact removal (0-1, higher better)
  - Measures how effectively batch effects are eliminated
  - Based on reduction in batch-associated variance (normalized)
  - 1.0 = perfect removal, 0.0 = artifacts remain

- **Stable↑**: Stability (0-1, higher better)
  - Measures consistency across random data subsets
  - Based on mean absolute error of repeated corrections (normalized)
  - 1.0 = perfectly stable, 0.0 = highly variable

- **DE↑**: Differential Expression preservation (-1 to 1, higher better)
  - Measures how well TB biomarker effect sizes are preserved
  - Based on correlation of Cohen's d values before/after correction
  - +1.0 = perfect preservation, 0.0 = no correlation, -1.0 = inverted

- **HK↑**: Housekeeping gene stability (0-1, higher better)
  - Measures how well stable genes remain stable across batches
  - Based on within-batch vs between-batch variance ratio
  - 1.0 = perfect stability, 0.0 = high batch variance

- **LinSep↑**: Linear separability (0-1, higher better)
  - Measures TB vs control classification performance (AUC)
  - Based on logistic regression with feature selection
  - 1.0 = perfect separation, 0.5 = random, 0.0 = inverted

- **Semantic↑**: Semantic anchoring quality (0-1, higher better)
  - Measures pathway-based sample matching accuracy
  - Based on whether pathway-similar samples have same TB status
  - 1.0 = perfect biological matching, 0.0 = random matching

- **Score↑**: Overall weighted score (0-1, higher better)
  - Weighted: 20% signal + 10% artifact + 5% stability + 25% DE + 15% HK + 15% LinSep + 10% Semantic
  - Signal and DE converted to 0-1 scale for scoring: (value+1)/2
  - Emphasizes biological signal preservation and classification performance

## Complete Analysis Output

```

📊 EXPERIMENT SUMMARY:
   Dataset: 10695 genes, 361 samples
   Methods tested: 12 configurations
   Validation: Signal preservation + 4 biological tests

🏆 BEST PERFORMER: posse_v50_balanced
   Version: POSSE v5.0 (ComBat-Initialized)
   Overall Score: 0.598
   Population Signal: 0.351, Individual Signal: 1.000
   Cross-Study Gen: 0.324, Classification: 0.947

📈 SIGNAL PRESERVATION INSIGHTS:
   • Population-level leaders: posse_v50_aggressive, posse_v50_balanced, posse_v50_conservative
   • Individual-level leaders: gmm_baseline, pace_v31_pure_local, posse_v50_conservative
   • No methods excel at both population and individual signal preservation

📊 METHOD-SPECIFIC INSIGHTS:
   • POSSE methods: Best overall score 0.598
     - Population preservation: 0.351
     - Individual preservation: 1.000
   • ComBat: Excellent HK stability (nan)
     - Population preservation: 0.162
     - Individual preservation: 0.999
   • GMM: Individual vs Population trade-off revealed
     - Population preservation: 0.009 (POOR - destroys compositional signal)
     - Individual preservation: 1.000 (maintains sample rankings)
     - This explains Snakemake success vs test script failure!
   • PACE methods: Balanced approach, best score 0.564
     - Population preservation: 0.309
     - Individual preservation: 1.000

🔄 CROSS-STUDY GENERALIZATION:
   • Top performers: v34_housekeeping_focused (0.428), v35_affine_focused (0.428), v22_aggressive_centered (0.411)

🔍 KEY DISCOVERY - GMM BEHAVIOR EXPLAINED:
   • GMM preserves individual sample relationships (good for classification)
   • GMM destroys population-level statistics (bad for meta-analysis)
   • This dual behavior explains the Snakemake vs test script discrepancy
```
