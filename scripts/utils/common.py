# Common utility functions moved from scripts/metrics/util.py

import os
import pandas as pd


def split_metadata_genes(df):
    # Combat-seq returns genes as integers, so splitting by type is not enough.
    # Other methods return genes as floats, so we need to split by metadata.
    # "meta_" has been prepended to all metadata columns, so we can split by that.
    metadata_cols = [col for col in df.columns if col.startswith("meta_")]
    genes = df.drop(columns=metadata_cols)
    metadata = df[metadata_cols]
    return metadata, genes


def split_into_batches(df, batch_col):
    _, genes = split_metadata_genes(df)
    batches = set(df[batch_col])
    return tuple((genes[df[batch_col] == batch] for batch in batches))