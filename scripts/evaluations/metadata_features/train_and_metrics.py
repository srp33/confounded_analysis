# Train one model per target and compute classification metrics

import os

# Set BLAS limits
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"

import re
import pandas as pd
import numpy as np
import argparse

from joblib import Parallel, delayed

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
# Parallelize across targets
# -------------------------


# -------------------------
# Classifier function
# -------------------------
def run_classifier(X_train, y_train, X_test, y_test, random_state=42):
    """Train classifier and compute metrics."""
    clf = HistGradientBoostingClassifier(max_iter=50, random_state=random_state)
    clf.fit(X_train, y_train)
    
    y_pred = clf.predict(X_test)
    y_proba = clf.predict_proba(X_test)

    n_classes = min(y_train.nunique(), y_test.nunique())
    
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
# Main
# -------------------------
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--csv", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--random_state", type=int, default=42)
    args = parser.parse_args()

    df = pd.read_csv(args.csv)
    adjuster, n_studies, test_source = parse_filename(args.csv)

    train_df = df[df["meta_source"].str.lower() != test_source.lower()]
    test_df = df[df["meta_source"].str.lower() == test_source.lower()]

    if train_df.empty or test_df.empty:
        raise ValueError("Train or test set is empty")

    feature_cols = [c for c in df.columns if not c.startswith("meta_")]

    target_cols = [
        "meta_er_status",
        "meta_her2_status",
        "meta_sex",
        "meta_chemotherapy",
        "meta_age_at_diagnosis_combined_lt50",
        "meta_age_at_diagnosis_combined_50_69",
        "meta_age_at_diagnosis_combined_ge70",
        "meta_menopause_status",
        "meta_histological_type",
    ]

    X_train_all = train_df[feature_cols]
    X_test_all = test_df[feature_cols]

    models_dir = os.path.join(args.outdir, adjuster, "models")
    os.makedirs(models_dir, exist_ok=True)

    def train_one_target(target):
        print("Training target: ", target)
        y_train = train_df[target]
        y_test = test_df[target]

        mask_train = y_train.notna()
        mask_test = y_test.notna()

        X_train = X_train_all.loc[mask_train]
        X_test = X_test_all.loc[mask_test]
        y_train = y_train.loc[mask_train]
        y_test = y_test.loc[mask_test]

        print(f"    Classes present in train dataset: {y_train.unique()}")
        print(f"    Classes present in test dataset: {y_test.unique()}")

        if y_train.nunique() < 2 or y_test.nunique() < 2:
            print(f"Skipping {target}: train or test has <2 classes")
            return None

        clf, metrics = run_classifier(
            X_train, y_train, X_test, y_test, args.random_state
        )

        try:
            import joblib
            joblib.dump(clf, os.path.join(models_dir, f"{target}.joblib"))
        except ImportError:
            pass

        return {
            "adjuster": adjuster,
            "n_studies": n_studies,
            "test_source": test_source,
            "target": target,
            **metrics,
        }

    n_jobs = int(os.environ.get("SLURM_CPUS_PER_TASK", 1))

    results = Parallel(n_jobs=n_jobs, backend="loky")(
        delayed(train_one_target)(t) for t in target_cols
    )

    metrics_df = pd.DataFrame([r for r in results if r is not None])

    out_dir = os.path.join(args.outdir, adjuster)
    os.makedirs(out_dir, exist_ok=True)
    metrics_df.to_csv(
        os.path.join(out_dir, f"{adjuster}_{n_studies}studies_test_{test_source}_metrics.csv"),
        index=False,
    )

    print("Done")

if __name__ == "__main__":
    main()