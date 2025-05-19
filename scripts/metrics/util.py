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

    def get_dataframe(self, path):
        try:
            return self.dataframes[path]
        except KeyError:
            self.dataframes[path] = pd.read_csv(path)
            return self.dataframes[path]

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

def cross_validate(df, predict_column, learner, iterations, folds, n_jobs, scale_numerics=False):
    meta_cols = list(df.select_dtypes(include=['object', 'int']).columns)

    X = df.drop(meta_cols, axis="columns")
    y = df[predict_column]

    if scale_numerics == True:
        X = robust_scale(X)

    #classes = []
    #for element in y:
    #    if (element not in classes) :
    #        classes.append(element)
    #y = label_binarize(y, classes = classes)

    scoring_metric = "roc_auc"

    scores = []
    for i in range(iterations):
        fit_params = learner[1]

        if "random_state" in fit_params:
            fit_params["random_state"] = i

        estimator = learner[0](**fit_params)

        kfold = StratifiedKFold(n_splits=folds, shuffle=True, random_state=i)
        iter_scores = list(cross_val_score(estimator, X, y, scoring=scoring_metric, cv=kfold, n_jobs=n_jobs))

        scores.append(sum(iter_scores) / len(iter_scores))

    return scores
