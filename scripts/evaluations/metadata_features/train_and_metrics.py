# Train one model per target and compute classification metrics

import os
import re
import json
import pandas as pd
import numpy as np
import argparse

from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.metrics import accuracy_score, roc_auc_score, confusion_matrix, matthews_corrcoef

import functools

print = functools.partial(print, flush=True)

# -------------------------
# Filename parsing
# -------------------------
def parse_filename(filename):
    basename = os.path.basename(filename).replace(".csv", "")

    pattern = r"^(?P<adjuster>.+)_(?P<n>\d+)studies_test_(?P<test>.+)$"
    match = re.match(pattern, basename)

    if not match:
        raise ValueError(f"Cannot parse filename: {basename}")

    return (
        match.group("adjuster"),
        int(match.group("n")),
        match.group("test")
    )

# -------------------------
# Classifier function
# -------------------------
def run_classifier(X_train, y_train, X_test, y_test, random_state=42):
    """Train classifier and compute metrics."""
    clf = HistGradientBoostingClassifier(max_iter=50, random_state=random_state)
    clf.fit(X_train, y_train)
    
    y_pred = clf.predict(X_test)
    y_proba = clf.predict_proba(X_test)

    n_classes = len(np.unique(y_train))
    metrics = {
        "n_classes": n_classes,
        "Accuracy": accuracy_score(y_test, y_pred),
        "MCC": matthews_corrcoef(y_test, y_pred)
    }
    
    # Binary Metrics
    if n_classes == 2:
        try:
            auc = roc_auc_score(y_test, y_proba[:, 1])
        except ValueError:
            auc = np.nan

        tn, fp, fn, tp = confusion_matrix(
            y_test, y_pred, labels=[0,1]
        ).ravel()

        metrics.update({
            "target_type":"binary",
            "ROC_AUC":auc,
            "Sensitivity":tp / (tp + fn) if (tp + fn) > 0 else np.nan,
            "Specificity": tn / (tn + fp) if (tn + fp) > 0 else np.nan,
            "TP": tp,
            "TN": tn,
            "FP": fp,
            "FN": fn,
        })
    
    # Multiclass metrics
    else:
        try:
            auc = roc_auc_score(
                y_test,
                y_proba,
                multi_class="ovr",
                average="macro"
            )
        except ValueError:
            auc = np.nan

        metrics.update({
            "target_type": "multiclass",
            "ROC_AUC_macro": auc,
        })

    return clf, metrics
    

# -------------------------
# Main function
# -------------------------
def main():
    parser = argparse.ArgumentParser(description="Run bootstrapped classifier on adjusted dataset")
    parser.add_argument("--csv", required=True, help="Input adjusted CSV file")
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument("--random_state", type=int, default=42)

    args = parser.parse_args()
    
    # -------------------------
    # Load data
    # -------------------------
    df = pd.read_csv(args.csv)
    adjuster, n_studies, test_source = parse_filename(args.csv)
    print(f"Adjuster: {adjuster}, n_studies: {n_studies}, test_source: {test_source}")
    
    # -------------------------
    # Split train/test
    # -------------------------
    train_df = df[df['meta_source'].str.lower() != test_source.lower()]
    test_df = df[df['meta_source'].str.lower() == test_source.lower()]
    
    if train_df.empty or test_df.empty:
        raise ValueError("Train or test set is empty — check test_source matching")
    
    feature_cols = [c for c in df.columns if not c.startswith("meta_")]

    target_cols = [
        'meta_er_status', 
        'meta_her2_status', 
        'meta_sex', 
        'meta_chemotherapy', 
        'meta_age_at_diagnosis_combined_lt50', 
        'meta_age_at_diagnosis_combined_50_69', 
        'meta_age_at_diagnosis_combined_ge70', 
        'meta_menopause_status',
        'meta_histological_type']

    X_train_all = train_df[feature_cols]
    X_test_all = test_df[feature_cols]

    metrics_rows = []
    models_dir = os.path.join(args.outdir, adjuster, "models")
    os.makedirs(models_dir, exist_ok=True)

    # Loop over targets
    for target in target_cols:
        print(f"Training target: {target}")

        y_train = train_df[target]
        y_test = train_df[target]

        # Drop NaNs PER TARGET
        train_mask = y_train.notna()
        test_mask = y_train.notna()

        X_train = X_train_all.loc[train_mask]
        X_test = X_test_all.loc[test_mask]
        y_train = y_train.loc[train_mask]
        y_test = y_test.loc[test_mask]

        # Skip degenerate targets
        if y_train.nunique() < 2:
            print(f"Skipping {target}: <2 classes")
            continue
        
        clf, metrics = run_classifier(
            X_train,
            y_train,
            X_test,
            y_test,
            random_state=args.random_state
        )

        metrics_row = {
            "adjuster": adjuster,
            "n_studies": n_studies,
            "test_source": test_source,
            "target": target,
            **metrics
        }
        metrics_rows.append(metrics_row)

        # Optional: save model
        model_path = os.path.join(models_dir, f"{target}.joblib")
        try:
            import joblib
            joblib.dump(clf, model_path)
        except ImportError:
            pass
 

    # -------------------------
    # Save results
    # -------------------------
    out_dir = os.path.join(args.outdir, adjuster)
    os.makedirs(out_dir, exist_ok=True)
    
    metrics_df = pd.DataFrame(metrics_rows)
    metrics_path = os.path.join(out_dir, f"{adjuster}_metrics.csv")
    metrics_df.to_csv(metrics_path, index=False)

    print(f"Saved metrics to: {metrics_path}")
    
if __name__ == "__main__":
    main()