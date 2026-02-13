#!/usr/bin/env python3

import os
import re
import argparse
import pandas as pd
import numpy as np
import joblib

from sklearn.inspection import permutation_importance

# -------------------------
# Filename parsing
# -------------------------
def parse_filename(filename):
    basename = os.path.basename(filename).replace(".csv", "")
    pattern = r"^(?P<adjuster>.+)_(?P<n>\d+)studies_test_(?P<test>.+)$"
    match = re.match(pattern, basename)
    if not match:
        raise ValueError(f"Cannot parse filename: {basename}")
    return match.group("adjuster"), int(match.group("n")), match.group("test")

# -------------------------
# Main
# -------------------------
def main():
    parser = argparse.ArgumentParser(description="Compute permutation importance per target")
    parser.add_argument("--csv", required=True, help="Input CSV")
    parser.add_argument("--models_dir", required=True, help="Directory containing trained models.")
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument("--n_repeats", type=int, default=3)
    parser.add_argument("--n_jobs", type=int, default=1)
    parser.add_argument("--random_state", type=int, default=42)
    args = parser.parse_args()

    # -------------------------
    # Load data
    # -------------------------
    df = pd.read_csv(args.csv)
    adjuster, n_studies, test_source = parse_filename(args.csv)

    train_df = df[df["meta_source"].str.lower() != test_source.lower()]
    test_df = df[df["meta_source"].str.lower() == test_source.lower()]

    if test_df.empty:
        raise ValueError("Test set is empty")

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

    X_test_all = test_df[feature_cols]

    print("======================================")
    print(f"Adjuster: {adjuster}")
    print(f"CSV: {args.csv}")
    print(f"Total samples: {df.shape[0]}")
    print(f"Test samples: {test_df.shape[0]}")
    print(f"Number of features: {len(feature_cols)}")
    print("======================================")

    if not os.path.isdir(args.models_dir):
        raise FileNotFoundError(f"Models directory not found: {args.models_dir}")

    # -------------------------
    # Permutation importance
    # -------------------------
    importance_dict = {}

    for target in target_cols:
        print(f"\nPermutation importance for target: {target}")

        model_path = os.path.join(args.models_dir, f"{target}.joblib")
        if not os.path.exists(model_path):
            print(f"  Model not found for {target}, skipping")
            continue

        clf = joblib.load(model_path)

        y_test = test_df[target]
        test_mask = y_test.notna()

        X_test = X_test_all.loc[test_mask]
        y_test = y_test.loc[test_mask]

        # Skip degenerate targets
        if y_test.nunique() < 2:
            print(f"  Skipping {target}: <2 classes in test set")
            continue

        model_classes = clf.classes_
        test_classes = np.unique(y_test)
        if set(test_classes) != set(model_classes):
            print(f"  Skipping {target}: model has classes {model_classes}, test has {test_classes}")
            continue

        scoring = "roc_auc" if len(clf.classes_) == 2 else "roc_auc_ovr"
        print(f"  Model type: {type(clf)}, scoring={scoring}")
        print(f"  Computing permutation importance...")

        r = permutation_importance(
            clf,
            X_test,
            y_test,
            scoring=scoring, 
            n_repeats=args.n_repeats,
            n_jobs=args.n_jobs,
            random_state=args.random_state,
        )

        importance_dict[target] = r.importances_mean

    # -------------------------
    # Save importance table
    # -------------------------
    os.makedirs(args.outdir, exist_ok=True)
    out_dir = os.path.join(args.outdir, "permutation_importance", adjuster)
    os.makedirs(out_dir, exist_ok=True)

    importance_df = pd.DataFrame(importance_dict, index=feature_cols)
    out_path = os.path.join(out_dir, f"{n_studies}_{test_source}_permutation_importance.csv")
    importance_df.to_csv(out_path)
    print(f"\nSaved permutation importance to: {out_path}")

if __name__ == "__main__":
    main()