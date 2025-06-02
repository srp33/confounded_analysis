import os
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


def split_metadata_genes(df):
    # Combat-seq returns genes as integers, so splitting by type is not enough.
    # Other methods return genes as floats, so I need to split by metadata.
    # I prepended "meta_" to all metadata columns, so I can split by that.
    metadata_cols = [col for col in df.columns if col.startswith("meta_")]
    genes = df.drop(columns=metadata_cols)
    metadata = df[metadata_cols]
    return metadata, genes


def split_into_batches(df, batch_col):
    _, genes = split_metadata_genes(df)
    batches = set(df[batch_col])
    return tuple((genes[df[batch_col] == batch] for batch in batches))



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