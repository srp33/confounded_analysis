import argparse
import os
from sklearn.ensemble import RandomForestClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.neighbors import NeighborhoodComponentsAnalysis
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.metrics import mutual_info_score


import os.path
from os import path
from pathlib import Path

import sys
import time
import random
from util import *
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input-dir", type=Path, help="Input directory", required=True)
parser.add_argument("-o", "--output-path", type=Path, help="Path to output file", required=True)
parser.add_argument("-b", "--batch-col", help="Batch column", required=True)
parser.add_argument("-c", "--column", help="Prediction column", required=True)
parser.add_argument("-m", "--confusion-path", type=Path, help="Path to confusion matrix file", required=False)
parser.add_argument('--class0', nargs='+', required=False, help='Values to map to class 0')
parser.add_argument('--class1', nargs='+', required=False, help='Values to map to class 1')
args = parser.parse_args()

cache = DataFrameCache()
dataset = os.path.basename(args.input_dir)
df = cache.get_dataframe(args.input_dir / "unadjusted.csv")

if not os.path.exists(args.output_path):
    with open(args.output_path, "w") as output_file:
        output_file.write("metric,adjuster,dataset,column,value\n")

if args.class0 and args.class1:
    # Map classes and flag invalid entries
    valid_classes = set(args.class0 + args.class1)
    df[args.column] = df[args.column].astype(str)

    invalid = df[~df[args.column].isin(valid_classes)]
    if not invalid.empty:
        print(f"Unmapped classes found: {invalid[args.column].unique()}")
    
    df[args.column] = df[args.column].apply(
        lambda x: 0 if x in args.class0 else 1 if x in args.class1 else np.nan
    ).dropna()  # Remove rows with unmapped classes
    
results = []

metric = "mutual_info_score"

# Print the confusion matrix between the two classes and the two batches
batch_col = df[args.batch_col]
predict_column = df[args.column]

print(f"Mutual Information Score between {args.column} and batch in {dataset} dataset:")
mutual_info = mutual_info_score(batch_col, predict_column) / np.log(2) # Convert to bits
print(mutual_info)
unique_batches = sorted(batch_col.unique())
unique_classes = sorted(predict_column.unique())

def calculate_counts(df, batch_col, predict_column):
    classes = df[predict_column].unique()
    batches = df[batch_col].unique()
    counts = [
        [
            len(df[(df[batch_col] == batch) & (df[predict_column] == cls)])
            for batch in batches
        ]
        for cls in classes
    ]
    counts = pd.DataFrame(counts, columns=classes, index=batches)
    return counts

counts = calculate_counts(df, args.batch_col, args.column)
print(counts)
print()

def product_of_marginals(counts):
    """
    Calculate the product of marginals for the given counts DataFrame.
    """
    row_sums = counts.sum(axis=1)
    col_sums = counts.sum(axis=0)
    total_sum = counts.values.sum()
    
    product = np.outer(row_sums, col_sums) / total_sum
    return pd.DataFrame(product, index=counts.index, columns=counts.columns)

product_counts = product_of_marginals(counts)
print("Product of Marginals:")
print(product_counts)
print()

def side_by_side_markdown(df1, df2, num_tabs=3):
    """
    Return two DataFrames as markdown side by side for comparison.
    Assume the same number of rows in both DataFrames.
    """
    df1_lines = df1.to_markdown().split('\n')
    df2_lines = df2.to_markdown().split('\n')
    max_length = max(len(line) for line in df1_lines + df2_lines)
    new_lines = []
    for line1, line2 in zip(df1_lines, df2_lines):
        line1 = line1.ljust(max_length)
        line2 = line2.ljust(max_length)
        new_lines.append(f"{line1}{' ' * num_tabs}{line2}")
    return '\n'.join(new_lines)



if args.confusion_path:
    with open(args.confusion_path, 'a') as conf_file:
        conf_file.write(f"#### Counts for Dataset: {dataset}\n\n")
        conf_file.write(f"Batches (row): {args.batch_col}\t\tPredictive Var (column): {args.column}\n\n")
        conf_file.write(f"**Mutual Information Score:** {round(mutual_info, 5)}\n\n")
        conf_file.write("#### Counts:\t\t\t\t\tProduct of Marginals:\n\n")
        conf_file.write(side_by_side_markdown(counts, product_counts))
        conf_file.write("\n\n---\n\n")

line = [metric, "Any", dataset, args.column, str(mutual_info)]

with open(args.output_path, 'a') as output_file:
    output_file.write(",".join(line) + "\n")
