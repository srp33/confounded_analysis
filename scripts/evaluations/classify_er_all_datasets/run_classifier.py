import os
import sys
import pandas as pd
import numpy as np
import argparse
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.metrics import (accuracy_score, roc_auc_score, confusion_matrix, matthews_corrcoef)

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
    model = HistGradientBoostingClassifier(max_iter=100, random_state=random_state)

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
    parser.add_argument("--bootstrap", type=int, default=1, help="Number of bootstrap runs (default: 1 = no bootstrap)")

    args = parser.parse_args()

    csv_file = args.csv
    output_dir = args.outdir
    n_boot = args.bootstrap

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

    print(f"[DEBUG] X_train shape: {X_train.shape}, y_train length: {len(y_train)}")
    print(f"[DEBUG] X_test shape: {X_test.shape}, y_test length: {len(y_test)}")

    # Assuming X_train is a DataFrame and y_train is a Series
    mask = y_train.notna()  # True for rows where y_train is not NaN
    X_train_clean = X_train[mask]
    y_train_clean = y_train[mask]

    mask_test = y_test.notna()
    X_test_clean = X_test[mask_test]
    y_test_clean = y_test[mask_test]

    # Bootstrapping logic
    results = []
    rng = np.random.default_rng(seed=123)

    for b in range(n_boot):
        print(f"\n=== Bootstrap iteration {b+1}/{n_boot} ===")
        if n_boot == 1 or b == 0:
            Xb, yb = X_train_clean, y_train_clean
        else:
            indices = rng.integers(0, len(X_train_clean), len(X_train_clean))
            Xb = X_train_clean.iloc[indices]
            yb = y_train_clean.iloc[indices]

        metrics = run_classifier(Xb, yb, X_test_clean, y_test_clean, random_state=42 + b)

        metrics['adjuster'] = adjuster
        metrics['subset_file'] = os.path.basename(csv_file)
        metrics['test_source'] = test_source
        metrics['bootstrap_iter'] = b + 1

        results.append(metrics)

    # --- Use the provided output directory ---
    results_dir = os.path.join(output_dir, adjuster)
    os.makedirs(results_dir, exist_ok=True)

    result_file = os.path.join(
        results_dir,
        f"{os.path.basename(csv_file).replace('.csv', '')}_metrics.csv"
    )
    
    df_results = pd.DataFrame(results)
    df_results.to_csv(result_file, index=False)

    print(f"\nSaved classifier metrics (including bootstrap runs): ")
    print(result_file)

if __name__ == "__main__":
    main()
