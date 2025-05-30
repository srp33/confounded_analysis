import os
import glob
import numpy as np
import pandas as pd
from sklearn.model_selection import cross_val_score
from sklearn.preprocessing import robust_scale
from sklearn.preprocessing import label_binarize
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import GridSearchCV
from pathlib import Path
from joblib import Parallel, delayed

def __sort_by_length(strings):
    return sorted(strings, key=len)

def __csvs_in_folder(folder):
    return glob.glob("{}/*.csv".format(folder))

def __get_filenames(paths):
    return [os.path.split(path)[-1] for path in paths]

def __infer_adjustment_dependencies(filenames):
    # When my adjusters are run, the output paths look like "base_adjuster.csv"
    # and the unadjusted one looks like "base.csv". So the shortest csv name
    # that is also contained in other csv names is probably the unadjusted csv.
    filenames = __sort_by_length(filenames)
    adjusted = []
    while (not adjusted) and filenames:
        unadjusted = filenames.pop(0)
        no_extension, _ = os.path.splitext(unadjusted)
        adjusted = [filename for filename in filenames if no_extension in filename]
    return {unadjusted: adjusted}

def __prepend_folder(path_dict, folder):
    for unadjusted, adjusted in path_dict.items():
        return {
            os.path.join(folder, unadjusted): [
                os.path.join(folder, filename) for filename in adjusted
            ]
        }

def get_dataset_path_dict(folder):
    """Get a dictionary like {unadjusted_path: [adjusted_path1, ...]} given a folder.
    """
    csv_paths = __csvs_in_folder(folder)
    filenames = __get_filenames(csv_paths)
    dependencies = __infer_adjustment_dependencies(filenames)
    return __prepend_folder(dependencies, folder)

def log_scale(df):
    # Get rid of negative values
    df = df.where(df.min() < 0, df - df.min())
    return np.log(df + 1.0)


class Logger(object):
    def __init__(self, metric):
        self.metric = metric
        self.values = {
            "adjuster": [],
            "dataset": [],
            "value": [],
        }

    def log(self, adjuster, dataset, value):
        self.values["adjuster"].append(adjuster)
        self.values["dataset"].append(dataset)
        self.values["value"].append(value)

    def save(self, path):
        # Appends new values to existing file
        df = pd.DataFrame(self.values)
        if os.path.exists(path):
            df = pd.concat([df, pd.read_csv(path)])
        df.to_csv(path, index=False)

    def save_pivoted(self, path):
        df = pd.DataFrame(self.values).drop("metric", axis=1)
        pivoted = df.pivot(index='dataset', columns='adjuster', values='value')
        pivoted.to_csv(Path(path) / self.metric)


class DataFrameCache(object):
    def __init__(self):
        self.dataframes = {} # path: dataframe

    def get_dataframe(self, file_path):
        if file_path in self.dataframes:
            return self.dataframes[file_path]
        else:
            try:
                df = pd.read_csv(file_path)
                self.dataframes[file_path] = df
                return df
            except FileNotFoundError:
                print(f"Process {os.getpid()}: ERROR - File not found: {file_path}", flush=True)
                return pd.DataFrame()
            except Exception as e:
                print(f"Process {os.getpid()}: ERROR - Could not read file {file_path} due to: {e}", flush=True)
                return pd.DataFrame()

def no_extension(path):
    filename = os.path.split(path)[-1]
    return os.path.splitext(filename)[0]

def split_discrete_continuous(df):
    discrete_types = ["int", "object"]
    discrete = df.select_dtypes(include=discrete_types)
    continuous = df.select_dtypes(exclude=discrete_types)
    return discrete, continuous

def split_into_batches(df, batch_col):
    discrete, continuous = split_discrete_continuous(df)
    batches = set(df[batch_col])
    return tuple((continuous[df[batch_col] == batch] for batch in batches))



def _run_single_iteration(iteration_idx, X_processed, y_processed, learner_config, num_folds, metric):
    # Make a copy of params, so we can modify random_state locally
    fit_params = learner_config.get("fit_params", {}).copy()
    if "random_state" in fit_params:
        fit_params["random_state"] = iteration_idx

    estimator = learner_config["algorithm"](**fit_params)
    kfold = StratifiedKFold(n_splits=num_folds, shuffle=True, random_state=iteration_idx)

    # n_jobs = 3 works here, currently running on 30 threads and averaging ~20 cpus of usage
    iter_scores = cross_val_score(estimator, X_processed, y_processed, scoring=metric, cv=kfold, n_jobs=3)

    if len(iter_scores) > 0:
        return sum(iter_scores) / len(iter_scores)
    return 0.0 # Return 0.0 if no scores are produced


def cross_validate(df, predict_column, learner, iterations, folds, n_jobs, scale_numerics=False):
    if df.empty:
        return []

    # Target variable
    y = df[predict_column].copy()

    # Remove target and other categorical columns
    X = df.drop(columns=[predict_column]).select_dtypes(exclude=['object', 'int']).copy()

    scoring_metric = "roc_auc"

    # Apply transformation if specified
    if "transform" in learner:
        transform_params = learner.get("transform_params", {})
        transformer = learner["transform"](**transform_params)
        transformer.fit(X, y) # Fit transformer
        X = transformer.transform(X) # Transform X

    # n_jobs for Parallel controls how many iterations run at once.
    scores = Parallel(n_jobs=n_jobs)(
        delayed(_run_single_iteration)(i, X, y, learner, folds, scoring_metric)
        for i in range(iterations)
    )

    return scores