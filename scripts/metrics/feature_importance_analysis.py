#!/usr/bin/env python3
"""
Feature Importance Analysis for HistGradientBoostingClassifier

This script analyzes feature importances of a HistGradientBoostingClassifier
trained on any dataset. It trains models on:
1. Each meta_source individually
2. Combined data from both sources
3. Cross-dataset evaluation (train on one source, test on another)

For feature importances, it uses:
1. sklearn.inspection.permutation_importance

Results are saved to CSV files for further analysis.
"""

import argparse
import os
import numpy as np
import pandas as pd
from pathlib import Path
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.inspection import permutation_importance
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_auc_score, accuracy_score, matthews_corrcoef, make_scorer
import matplotlib.pyplot as plt
import seaborn as sns
import time
import warnings
warnings.filterwarnings('ignore')

def print_now(*args, **kwargs):
    """Print with immediate flush"""
    print(*args, flush=True, **kwargs)

def load_and_prepare_data(data_path):
    """Load data and prepare features and target"""
    print_now(f"Loading data from {data_path}")
    df = pd.read_csv(data_path)
    
    # Separate features and target
    target_col = 'meta_er_status'
    source_col = 'meta_source'
    
    if target_col not in df.columns:
        raise ValueError(f"Target column '{target_col}' not found in data")
    if source_col not in df.columns:
        raise ValueError(f"Source column '{source_col}' not found in data")
    
    # Get feature columns (exclude meta columns)
    feature_cols = [col for col in df.columns if not col.startswith('meta_')]
    X = df[feature_cols].select_dtypes(include=[np.number])
    y = df[target_col]
    sources = df[source_col]
    
    print_now(f"Data shape: {df.shape}")
    print_now(f"Features: {X.shape[1]} numeric columns")
    print_now(f"Target distribution: {y.value_counts().to_dict()}")
    print_now(f"Sources: {sources.unique()}")
    
    return X, y, sources, feature_cols

def train_and_analyze_model(X_train, y_train, X_test, y_test, model_name, random_state=42, fast_mode=False):
    """Train model and extract permutation feature importances"""
    print_now(f"Training {model_name}...")
    start_time = time.time()
    
    # Train model
    model = HistGradientBoostingClassifier(
        max_iter=100, 
        random_state=random_state,
        early_stopping=True,
        validation_fraction=0.1,
        n_iter_no_change=10
    )
    
    model.fit(X_train, y_train)
    train_time = time.time() - start_time
    print_now(f"  Training completed in {train_time:.2f} seconds")
    
    # Get predictions and performance
    y_pred = model.predict(X_test)
    y_proba = model.predict_proba(X_test)[:, 1]
    
    accuracy = accuracy_score(y_test, y_pred)
    mcc = matthews_corrcoef(y_test, y_pred)
    try:
        auc = roc_auc_score(y_test, y_proba)
    except ValueError:
        auc = np.nan  # Handle case where only one class is present
    
    print_now(f"{model_name} - Accuracy: {accuracy:.4f}, MCC: {mcc:.4f}, AUC: {auc:.4f}")
    
    # Get permutation importances using MCC as the primary metric
    print_now(f"Computing permutation importances for {model_name}...")
    print_now(f"  Features: {X_test.shape[1]}, Test samples: {X_test.shape[0]}")
    
    perm_start_time = time.time()
    mcc_scorer = make_scorer(matthews_corrcoef)
    
    # Optimize parameters based on dataset size and fast_mode
    n_features = X_test.shape[1]
    n_samples = X_test.shape[0]
    
    
    perm_importance = permutation_importance(
        model, X_test, y_test, 
        n_repeats=3,
        random_state=random_state,
        n_jobs=-1,
        scoring=mcc_scorer
    )
    
    perm_time = time.time() - perm_start_time
    print_now(f"  Permutation importance completed in {perm_time:.2f} seconds")
    
    results = {
        'model_name': model_name,
        'accuracy': accuracy,
        'mcc': mcc,
        'auc': auc,
        'perm_importances_mean': perm_importance.importances_mean,
        'perm_importances_std': perm_importance.importances_std,
        'feature_names': X_train.columns.tolist(),
        'train_time': train_time,
        'perm_time': perm_time
    }
    
    return results

def analyze_single_source(X, y, sources, source_name, output_dir):
    """Analyze feature importances for a single source"""
    print_now(f"\n{'='*60}")
    print_now(f"ANALYZING SOURCE: {source_name}")
    print_now(f"{'='*60}")
    
    # Filter data for this source
    source_mask = sources == source_name
    X_source = X[source_mask]
    y_source = y[source_mask]
    
    print_now(f"Source {source_name} - Samples: {len(X_source)}, Target distribution: {y_source.value_counts().to_dict()}")
    
    # Check if we have enough samples and both classes
    if len(X_source) < 10:
        print_now(f"Warning: Too few samples ({len(X_source)}) for source {source_name}. Skipping.")
        return None
    
    if len(y_source.unique()) < 2:
        print_now(f"Warning: Only one class present in source {source_name}. Skipping.")
        return None
    
    # Split data
    X_train, X_test, y_train, y_test = train_test_split(
        X_source, y_source, 
        test_size=0.3, 
        random_state=42, 
        stratify=y_source
    )
    
    # Train and analyze
    results = train_and_analyze_model(
        X_train, y_train, X_test, y_test, 
        f"Source: {source_name}", 
        random_state=42
    )
    
    # Save results
    save_importance_results(results, output_dir, f"source_{source_name}")
    
    return results

def analyze_combined_sources(X, y, sources, output_dir):
    """Analyze feature importances for combined data"""
    print_now(f"\n{'='*60}")
    print_now(f"ANALYZING COMBINED SOURCES")
    print_now(f"{'='*60}")
    
    print_now(f"Combined data - Samples: {len(X)}, Target distribution: {y.value_counts().to_dict()}")
    
    # Split data
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, 
        test_size=0.3, 
        random_state=42, 
        stratify=y
    )
    
    # Train and analyze
    results = train_and_analyze_model(
        X_train, y_train, X_test, y_test, 
        "Combined Sources", 
        random_state=42
    )
    
    # Save results
    save_importance_results(results, output_dir, "combined_sources")
    
    return results

def analyze_cross_dataset(X, y, sources, source1, source2, output_dir):
    """Analyze feature importances for cross-dataset evaluation (train on source1, test on source2)"""
    print_now(f"\n{'='*60}")
    print_now(f"CROSS-DATASET ANALYSIS: Train on {source1}, Test on {source2}")
    print_now(f"{'='*60}")
    
    # Get data for each source
    source1_mask = sources == source1
    source2_mask = sources == source2
    
    X_train = X[source1_mask]
    y_train = y[source1_mask]
    X_test = X[source2_mask]
    y_test = y[source2_mask]
    
    print_now(f"Train ({source1}) - Samples: {len(X_train)}, Target distribution: {y_train.value_counts().to_dict()}")
    print_now(f"Test ({source2}) - Samples: {len(X_test)}, Target distribution: {y_test.value_counts().to_dict()}")
    
    # Check if we have enough samples and both classes in both datasets
    if len(X_train) < 10 or len(X_test) < 10:
        print_now(f"Warning: Too few samples for cross-dataset analysis. Skipping.")
        return None
    
    if len(y_train.unique()) < 2 or len(y_test.unique()) < 2:
        print_now(f"Warning: Only one class present in one of the datasets. Skipping.")
        return None
    
    # Train and analyze
    results = train_and_analyze_model(
        X_train, y_train, X_test, y_test, 
        f"Train on {source1}, Test on {source2}", 
        random_state=42
    )
    
    # Save results
    save_importance_results(results, output_dir, f"cross_train_{source1}_test_{source2}")
    
    return results

def save_importance_results(results, output_dir, prefix):
    """Save feature importance results to CSV files"""
    feature_names = results['feature_names']
    
    # Create DataFrame with importance information
    importance_df = pd.DataFrame({
        'feature': feature_names,
        'perm_importance_mean': results['perm_importances_mean'],
        'perm_importance_std': results['perm_importances_std']
    })
    
    # Add model metadata
    importance_df['model_name'] = results['model_name']
    importance_df['model_accuracy'] = results['accuracy']
    importance_df['model_mcc'] = results['mcc']
    importance_df['model_auc'] = results['auc']
    
    # Sort by permutation importance
    importance_df = importance_df.sort_values('perm_importance_mean', ascending=False)
    
    # Save to CSV
    csv_path = output_dir / f"{prefix}_feature_importances.csv"
    importance_df.to_csv(csv_path, index=False, float_format='%.6f')
    print_now(f"Saved feature importances to: {csv_path}")
    
    # Print top 10 features
    print_now(f"\nTop 10 features by permutation importance (MCC-based) for {results['model_name']}:")
    top_features = importance_df.head(10)
    for idx, row in top_features.iterrows():
        print_now(f"  {row['feature']}: perm={row['perm_importance_mean']:.6f}±{row['perm_importance_std']:.6f}")

def create_comparison_plots(all_results, output_dir):
    """Create comparison plots of feature importances across different models"""
    print_now(f"\n{'='*60}")
    print_now("CREATING COMPARISON PLOTS")
    print_now(f"{'='*60}")
    
    if not all_results:
        print_now("No results to plot")
        return
    
    # Combine all results into a single DataFrame
    combined_data = []
    for result in all_results:
        if result is None:
            continue
        
        for i, feature in enumerate(result['feature_names']):
            combined_data.append({
                'model': result['model_name'],
                'feature': feature,
                'perm_importance_mean': result['perm_importances_mean'][i],
                'perm_importance_std': result['perm_importances_std'][i]
            })
    
    if not combined_data:
        print_now("No data to plot")
        return
    
    df_plot = pd.DataFrame(combined_data)
    
    # Get top 20 features by average permutation importance across all models
    top_features = (df_plot.groupby('feature')['perm_importance_mean']
                   .mean()
                   .sort_values(ascending=False)
                   .head(20)
                   .index.tolist())
    
    df_plot_top = df_plot[df_plot['feature'].isin(top_features)]
    
    # Create comparison plots
    fig, axes = plt.subplots(1, 2, figsize=(16, 8))
    
    # Plot 1: Permutation importance comparison heatmap
    pivot_perm = df_plot_top.pivot(index='feature', columns='model', values='perm_importance_mean')
    pivot_perm = pivot_perm.reindex(top_features)  # Maintain order
    
    sns.heatmap(pivot_perm, annot=True, fmt='.4f', cmap='viridis', ax=axes[0])
    axes[0].set_title('Permutation Importance Comparison - MCC-based (Top 20 Features)')
    axes[0].set_xlabel('Model')
    axes[0].set_ylabel('Feature')
    
    # Plot 2: Top features bar plot
    top_10_features = top_features[:10]
    df_top_10 = df_plot[df_plot['feature'].isin(top_10_features)]
    
    # Create grouped bar plot
    x_pos = np.arange(len(top_10_features))
    width = 0.8 / len(df_plot['model'].unique())
    models = df_plot['model'].unique()
    
    for i, model in enumerate(models):
        model_data = df_top_10[df_top_10['model'] == model]
        model_data = model_data.set_index('feature').reindex(top_10_features)
        
        axes[1].bar(x_pos + i*width, 
                   model_data['perm_importance_mean'], 
                   width, 
                   label=model,
                   yerr=model_data['perm_importance_std'],
                   capsize=3)
    
    axes[1].set_xlabel('Feature')
    axes[1].set_ylabel('Permutation Importance')
    axes[1].set_title('Top 10 Features - Permutation Importance by Model')
    axes[1].set_xticks(x_pos + width * (len(models) - 1) / 2)
    axes[1].set_xticklabels(top_10_features, rotation=45, ha='right')
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    # Save plot
    plot_path = output_dir / "feature_importance_comparison.png"
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    plt.close()
    print_now(f"Saved comparison plot to: {plot_path}")

def main():
    parser = argparse.ArgumentParser(description="Analyze feature importances of HistGradientBoostingClassifier")
    parser.add_argument('--data-path', 
                       help='Path to the input data file')
    parser.add_argument('--output-dir', 
                       default='/outputs/metrics/feature_importance',
                       help='Output directory for results')
    parser.add_argument('--create-plots', 
                       action='store_true',
                       help='Create comparison plots')
    
    args = parser.parse_args()
    
    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print_now("="*80)
    print_now("FEATURE IMPORTANCE ANALYSIS FOR HISTGRADIENTBOOSTINGCLASSIFIER")
    print_now("="*80)
    print_now(f"Data path: {args.data_path}")
    print_now(f"Output directory: {output_dir}")
    
    # Load and prepare data
    try:
        X, y, sources, feature_cols = load_and_prepare_data(args.data_path)
    except Exception as e:
        print_now(f"Error loading data: {e}")
        return
    
    # Store all results for comparison
    all_results = []
    
    # Analyze each source individually
    unique_sources = sources.unique()
    for source in unique_sources:
        result = analyze_single_source(X, y, sources, source, output_dir)
        if result:
            all_results.append(result)
    
    # Analyze combined sources
    combined_result = analyze_combined_sources(X, y, sources, output_dir)
    if combined_result:
        all_results.append(combined_result)
    
    # Analyze cross-dataset scenarios (train on one source, test on another)
    if len(unique_sources) >= 2:
        for i, source1 in enumerate(unique_sources):
            for j, source2 in enumerate(unique_sources):
                if i != j:  # Don't train and test on the same source
                    cross_result = analyze_cross_dataset(X, y, sources, source1, source2, output_dir)
                    if cross_result:
                        all_results.append(cross_result)
    
    # Create summary comparison
    if all_results:
        print_now(f"\n{'='*60}")
        print_now("SUMMARY COMPARISON")
        print_now(f"{'='*60}")
        
        for result in all_results:
            print_now(f"{result['model_name']}: Accuracy={result['accuracy']:.4f}, MCC={result['mcc']:.4f}, AUC={result['auc']:.4f}")
        
        # Save combined results
        combined_csv_path = output_dir / "all_feature_importances_combined.csv"
        combined_data = []
        
        for result in all_results:
            for i, feature in enumerate(result['feature_names']):
                combined_data.append({
                    'model': result['model_name'],
                    'feature': feature,
                    'perm_importance_mean': result['perm_importances_mean'][i],
                    'perm_importance_std': result['perm_importances_std'][i],
                    'model_accuracy': result['accuracy'],
                    'model_mcc': result['mcc'],
                    'model_auc': result['auc']
                })
        
        combined_df = pd.DataFrame(combined_data)
        combined_df.to_csv(combined_csv_path, index=False, float_format='%.6f')
        print_now(f"Saved combined results to: {combined_csv_path}")
        
        # Create plots if requested
        if args.create_plots:
            create_comparison_plots(all_results, output_dir)
    
    print_now(f"\n{'='*60}")
    print_now("ANALYSIS COMPLETE")
    print_now(f"{'='*60}")
    print_now(f"Results saved to: {output_dir}")

if __name__ == "__main__":
    main()