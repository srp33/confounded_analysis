import os
import glob
import numpy as np
import pandas as pd

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
            "metric": [],
            "adjuster": [],
            "dataset": [],
            "value": [],
        }

    def log(self, adjuster, dataset, value):
        self.values["metric"].append(self.metric)
        self.values["adjuster"].append(adjuster)
        self.values["dataset"].append(dataset)
        self.values["value"].append(value)

    def save(self, path):
        df = pd.DataFrame(self.values)
        df.to_csv(path, index=False)

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
