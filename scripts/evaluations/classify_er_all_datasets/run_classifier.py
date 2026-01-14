import os
import sys
import pandas as pd
import numpy as np
import argparse
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.metrics import (accuracy_score, roc_auc_score, confusion_matrix, matthews_corrcoef)
from sklearn.utils import resample
from sklearn.feature_selection import VarianceThreshold

import functools
print = functools.partial(print, flush=True)

def parse_test_source(filename):
    """Extract the test source from filename."""
    basename = os.path.basename(filename)
    parts = basename.split("_test_")
    if len(parts) != 2:
        raise ValueError(f"Cannot parse test source from filename: {basename}")
    test_source = parts[1].replace(".csv", "")
    return test_source

def run_classifier(X_train, y_train, X_test, y_test, random_state=42):
    model = HistGradientBoostingClassifier(max_iter=50, random_state=random_state)

    # Now fit the model
    model.fit(X_train, y_train)

    y_pred = model.predict(X_test)
    y_proba = model.predict_proba(X_test)[:, 1]

    tn, fp, fn, tp = confusion_matrix(y_test, y_pred, labels=[0,1]).ravel()
    acc = accuracy_score(y_test, y_pred)
    try: 
        auc = roc_auc_score(y_test, y_proba)
    except ValueError: 
        auc = np.nan
    sens = tp / (tp + fn) if (tp + fn) > 0 else np.nan
    spec = tn / (tn + fp) if (tn + fp) > 0 else np.nan
    mcc = matthews_corrcoef(y_test, y_pred)

    return {
        "Accuracy": acc,
        "ROC_AUC": auc,
        "Sensitivity": sens,
        "Specificity": spec,
        "MCC": mcc,
        "True Negative": tn,
        "False Positive": fp,
        "False Negative": fn,
        "True Positive": tp
    }

def main():
    parser = argparse.ArgumentParser(
        description="Train/test classifier with bootstrapping."
    )

    parser.add_argument("--csv", required=True, help="Input adjusted CSV file")
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument("--n_hvg", type=int, default=1000,
                    help="Number of highly variable genes to select")
    parser.add_argument("--chunk", type=int, default=0, help="Bootstrap chunk index (0-based)")
    parser.add_argument("--chunk-size", type=int, default=10, help="How many bootstraps per job")
    
    args = parser.parse_args()

    csv_file = args.csv
    output_dir = args.outdir
    n_hvg = args.n_hvg
    chunk = args.chunk
    chunk_size = args.chunk_size

    if not os.path.exists(csv_file):
        raise FileNotFoundError(f"CSV not found: {csv_file}")

    df = pd.read_csv(csv_file)

    if 'meta_er_status' not in df.columns or 'meta_source' not in df.columns:
        raise ValueError("CSV must contain 'meta_er_status' and 'meta_source' columns")

    # print(f"[DEBUG] Missing values in 'meta_er_status': {df['meta_er_status'].isna().sum()} / {len(df)} total rows")
    # print(f"[DEBUG] Unique meta_source values: {df['meta_source'].unique()[:10]}")
    # print(f"[DEBUG] Unique meta_er_status values: {df['meta_er_status'].unique()[:10]}")

    test_source = parse_test_source(csv_file)
    adjuster = os.path.basename(os.path.dirname(csv_file))

    # ✅ Case-insensitive split
    train_df = df[df['meta_source'].str.lower() != test_source.lower()]
    test_df = df[df['meta_source'].str.lower() == test_source.lower()]
 
    if train_df.empty or test_df.empty:
        print("[ERROR] Train or test set is empty — possible mismatch between test_source and meta_source values")
        print(f"[DEBUG] First few meta_source values: {df['meta_source'].head()}")
        sys.exit(1)

    # # ✅ Print unique label values to catch string vs numeric issues
    # print("[DEBUG] Unique values in meta_er_status (train):", train_df['meta_er_status'].unique())
    # print("[DEBUG] Unique values in meta_er_status (test):", test_df['meta_er_status'].unique())

    # # ✅ Print all meta columns for reference
    # meta_cols = [c for c in df.columns if c.startswith("meta_")]
    # print(f"[DEBUG] Meta columns found: {meta_cols}")
    # print("[DEBUG] Example meta data (first 5 rows):")
    # print(df[meta_cols].head())

    feature_cols = [c for c in df.columns if not c.startswith("meta_")]
    X_train, y_train = train_df[feature_cols], train_df['meta_er_status']
    X_test, y_test = test_df[feature_cols], test_df['meta_er_status']

    # Select top 1000 highly variable genes
    gene_variances = X_train.var(axis=0)

    top_genes = gene_variances.sort_values(ascending=False).head(n_hvg).index

    print(f"[INFO] Selecting top {len(top_genes)} highly variable genes")

    X_train = X_train[top_genes]
    X_test = X_test[top_genes]

    # Assuming X_train is a DataFrame and y_train is a Series
    mask = y_train.notna()  # True for rows where y_train is not NaN
    X_train_clean = X_train[mask]
    y_train_clean = y_train[mask]

    mask_test = y_test.notna()
    X_test_clean = X_test[mask_test]
    y_test_clean = y_test[mask_test]

    # Bootstrapping logic
    results = []
    chunk_seed_base = 10_000 + chunk * chunk_size
    rng = np.random.default_rng(chunk_seed_base)

    for local_iter in range(chunk_size):
        global_iter = chunk * chunk_size + local_iter
        bootstrap_seed = int(rng.integers(1_000_000_000)) # seed for this iteration

        print(f"\n=== Bootstrap global iter {global_iter} (chunk {chunk}, local {local_iter}) ===")

        X_boot, y_boot = resample(
            X_train_clean,
            y_train_clean,
            replace=True,
            n_samples=len(X_train_clean),
            stratify=y_train_clean,        # maintains class balance
            random_state=bootstrap_seed
        )

        metrics = run_classifier(
            X_boot,
            y_boot,
            X_test_clean,
            y_test_clean,
            random_state=bootstrap_seed
        )

        metrics["bootstrap_global"] = global_iter
        metrics["bootstrap_chunk"] = chunk
        metrics["bootstrap_local"] = local_iter
        metrics["bootstrap_seed"] = bootstrap_seed
        metrics["adjuster"] = adjuster
        metrics["subset_file"] = os.path.basename(csv_file)
        metrics["test_source"] = test_source

        results.append(metrics)

    # --- Use the provided output directory ---
    results_dir = os.path.join(output_dir, adjuster)
    os.makedirs(results_dir, exist_ok=True)

    result_file = os.path.join(
        results_dir,
        f"{os.path.basename(csv_file).replace('.csv', '')}_chunk{chunk}_metrics.csv"
    )
    
    df_results = pd.DataFrame(results)
    df_results.to_csv(result_file, index=False)

    print(f"\nSaved classifier metrics (including bootstrap runs): ")
    print(result_file)

if __name__ == "__main__":
    main()
