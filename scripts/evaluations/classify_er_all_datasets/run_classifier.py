import os
import sys
import pandas as pd
import numpy as np
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.metrics import (accuracy_score, roc_auc_score, confusion_matrix, matthews_corrcoef)

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
        "ROC AUC": auc,
        "Sensitivity": sens,
        "Specificity": spec,
        "MCC": mcc,
        "True Negative": tn,
        "False Positive": fp,
        "False Negative": fn,
        "True Positive": tp
    }

def main():
    if len(sys.argv) != 2:
        print("Usage: python run_classifier.py <adjusted_csv>")
        sys.exit(1)

    csv_file = sys.argv[1]
    if not os.path.exists(csv_file):
        print(f"File not found: {csv_file}")
        sys.exit(1)

    df = pd.read_csv(csv_file)
    test_source = parse_test_source(csv_file)
    adjuster = os.path.basename(os.path.dirname(csv_file))

    if 'meta_er_status' not in df.columns or 'meta_source' not in df.columns:
        raise ValueError("CSV must contain 'meta_er_status' and 'meta_source' columns")
    
    # Train/test split
    train_df = df[df['meta_source'] != test_source]
    test_df = df[df['meta_source'] == test_source]
    if train_df.empty or test_df.empty:
        raise ValueError(f"No train/test samples for test source {test_source}")
    
    # Features
    feature_cols = [c for c in df.columns if not c.startswith("meta_")]
    X_train, y_train = train_df[feature_cols].values, train_df['meta_er_status'].values
    X_test, y_test = test_df[feature_cols].values, test_df['meta_er_status'].values

    metrics = run_classifier(X_train, y_train, X_test, y_test)
    metrics['adjuster'] = adjuster
    metrics['subset_file'] = os.path.basename(csv_file)
    metrics['test_source'] = test_source

    # Save temporary per-job CSV
    results_dir = os.path.join("results", adjuster)
    os.makedirs(results_dir, exist_ok=True)
    result_file = os.path.join(results_dir, f"{os.path.basename(csv_file).replace('.csv','')}_metrics.csv")
    pd.DataFrame([metrics]).to_csv(result_file, index=False)

    print(f"Saved classifier metrics: {result_file}")

if __name__ == "__main__":
    main()
