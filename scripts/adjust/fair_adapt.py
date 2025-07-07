import numpy as np
import pandas as pd
from argparse import ArgumentParser
from aif360.sklearn.preprocessing import FairAdapt

# Parse command line args --------------------------
parser = ArgumentParser(description="AutoClass Imputation Example")
parser.add_argument("-i", "--input-file", help="Path to input CSV file", required=True)
parser.add_argument("-o", "--output-file", help="Path to output CSV file", required=True)
parser.add_argument("-b", "--batch-column", help="Column name for batch information", required=True)
args = parser.parse_args()

# Load the dataset --------------------------------
input_file = args.input_file
if not input_file.endswith('.csv'):
    raise ValueError("Input file must be a CSV file.")

try:
    df = pd.read_csv(input_file)
except Exception as e:
    raise ValueError(f"AutoClass Error reading the input file: {input_file}\n{e}")

# Perform the adjustments --------------------------

# Check for negative values. Assuming that microarray data is mean 0, while RNA-seq data is counts
# Make sure to filter out the metadata columns if they exist. Metadata columns start with 'meta_'
meta_cols = [col for col in df.columns if col.startswith('meta_')]
print(f"Found metadata columms: {meta_cols}")
gene_cols = [col for col in df.columns if col not in meta_cols]

batches = df[args.batch_column].values


# construct an adjacency matrix
adj_mat = pd.DataFrame(
    np.zeros((len(df.columns), len(df.columns)), dtype=int),
    index=df.columns.values,
    columns=df.columns.values
)

# Construct the adjacency matrix of the causal graph
# The metavariables cause the gene variables
adj_mat.loc[meta_cols, gene_cols] = 1

# Batch is the protected attribute
FA = FairAdapt(prot_attr=args.batch_column, adj_mat=adj_mat)

x = df[gene_cols + [args.batch_column]]
meta_without_batch = [col for col in meta_cols if col != args.batch_column]
y = df[meta_without_batch]
# Get first row
first_row = x.iloc[0, :]
x, y, _ = FA.fit_transform(x, y, first_row)
df =  pd.DataFrame(x, columns=gene_cols + [args.batch_column])
# Add metadata columns back to the DataFrame
df[meta_without_batch] = y

# Save the imputed data to a CSV file --------------------
output_file = args.output_file
try:
    df.to_csv(output_file, index=False)
    print(f"Imputed data saved to {output_file}")
except Exception as e:
    raise ValueError(f"AutoClass Error saving the output file: {output_file}\n{e}")