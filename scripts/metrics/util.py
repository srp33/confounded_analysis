# util.py - Machine learning utilities

import os
import sys
import glob
import numpy as np
import pandas as pd
from sklearn.model_selection import cross_validate
from sklearn.preprocessing import robust_scale
from sklearn.preprocessing import label_binarize
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import GridSearchCV
from pathlib import Path
from joblib import Parallel, delayed


def _run_single_iteration(iteration_idx, X, y, learner_config, n_folds, metrics):
    # Make a copy of params, so we can modify random_state locally
    fit_params = learner_config.get("fit_params", {}).copy()
    if "random_state" in fit_params:
        fit_params["random_state"] = iteration_idx

    estimator = learner_config["algorithm"](**fit_params)
    kfold = StratifiedKFold(n_splits=n_folds, shuffle=True, random_state=iteration_idx)

    scores = cross_validate(estimator, X, y, scoring=metrics, cv=kfold, n_jobs=3)

    test_names = {metric: "test" + "_" + metric for metric in metrics}
    average_scores = {metric: np.mean(scores[test]) for metric, test in test_names.items()}
    return average_scores


def repeated_cross_val(df, predict_column, learner, iterations, n_folds, n_jobs, metrics, scale_numerics=False):
    if df.empty:
        return []

    # Target variable
    if predict_column not in df.columns:
        print(f"Predict column {predict_column} not found in dataframe.")
        print(f"Columns beginning with 'meta_' are: {df.columns[df.columns.str.startswith('meta_')]}")
    y = df[predict_column].copy()

    # Remove target and other categorical columns
    X = df.drop(columns=[predict_column]).select_dtypes(exclude=['object', 'int']).copy()

    # Apply transformation if specified
    if "transform" in learner:
        transform_params = learner.get("transform_params", {})
        transformer = learner["transform"](**transform_params)
        transformer.fit(X, y) # Fit transformer
        X = transformer.transform(X) # Transform X

    # Parallel n_jobs controls how many iterations run at once.
    scores = Parallel(n_jobs=n_jobs)(
        delayed(_run_single_iteration)(i, X, y, learner, n_folds, metrics)
        for i in range(iterations)
    )

    return scores
