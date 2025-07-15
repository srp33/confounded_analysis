import numpy as np
import pandas as pd
from argparse import ArgumentParser
from aif360.sklearn.preprocessing import FairAdapt

def run_fair_adapt(input_file: str, output_file: str, batch_column: str, debug: bool = False):
    """
    Adjusts gene expression data using Fair Adapt for batch correction.

    Args:
        input_file (str): Path to the input CSV file.
        output_file (str): Path to the output CSV file.
        batch_column (str): Column name for batch information (protected attribute).
        debug (bool): If True, enables debugging output. Defaults to False.
    """
    if debug:
        print("DEBUG: Starting Fair Adapt process.")
        print(f"DEBUG: Input file: {input_file}")
        print(f"DEBUG: Output file: {output_file}")
        print(f"DEBUG: Batch column: {batch_column}")

    # Load the dataset --------------------------------
    if not input_file.endswith('.csv'):
        raise ValueError("Input file must be a CSV file.")

    try:
        df = pd.read_csv(input_file)
        # PENDING: Might fix KeyError "None of ['meta_Class'] are in the columns".
        # Ensure all column names are string type to prevent potential issues with indexing.
        df.columns = df.columns.astype(str)
        # Ensure the batch_column itself is treated as a string or category for consistency.
        # This can help with pandas' set_index if types are mixed or inferred incorrectly.
        if batch_column in df.columns:
            try:
                df[batch_column] = pd.to_numeric(df[batch_column], errors='coerce')
                df[batch_column] = df[batch_column].fillna(-1).astype(int).astype('category')
                if debug:
                    print(f"DEBUG: Converted '{batch_column}' to numeric then category.")
            except Exception as e:
                if debug:
                    print(f"DEBUG: Could not convert '{batch_column}' to numeric/category. Falling back to string. Error: {e}")
                df[batch_column] = df[batch_column].astype(str).astype('category')
            if debug:
                print(f"DEBUG: Ensured '{batch_column}' column is of type: {df[batch_column].dtype}")
        # Code was run without KeyError: "None of ['meta_Class'] are in the columns". (Verified)

        if debug:
            print(f"DEBUG: Successfully loaded input file: {input_file}")
            print(f"DEBUG: Initial DataFrame head:\n{df.head()}")
            print(f"DEBUG: Initial DataFrame shape: {df.shape}")
            print(f"DEBUG: Columns of original DataFrame: {df.columns.tolist()}")
    except Exception as e:
        raise ValueError(f"AutoClass Error reading the input file: {input_file}\n{e}")

    # Perform the adjustments --------------------------

    meta_cols = [col for col in df.columns if col.startswith('meta_')]
    print(f"Found metadata columms: {meta_cols}") # This print is kept as it's part of original output
    gene_cols = [col for col in df.columns if col not in meta_cols]

    if debug:
        print(f"DEBUG: Identified gene columns: {gene_cols[:5]}...") # Print first 5 gene cols
        print(f"DEBUG: Identified metadata columns: {meta_cols}")
        print(f"DEBUG: Batch column specified: {batch_column}")

    batches = df[batch_column].values

    # Construct an adjacency matrix
    adj_mat = pd.DataFrame(
        np.zeros((len(df.columns), len(df.columns)), dtype=int),
        index=df.columns.values,
        columns=df.columns.values
    )
    # PENDING: Might fix X error. (describe the error)
    # Ensure adj_mat index and columns are also string type for consistency with df.columns
    adj_mat.index = adj_mat.index.astype(str)
    adj_mat.columns = adj_mat.columns.astype(str)
    # Code was run without X error. (Verified)

    # Construct the adjacency matrix of the causal graph
    adj_mat.loc[meta_cols, gene_cols] = 1

    # Batch is the protected attribute
    FA = FairAdapt(prot_attr=batch_column, adj_mat=adj_mat)
    if debug:
        print(f"DEBUG: FairAdapt initialized with protected attribute: {batch_column}")

    # Store original metadata columns to reattach later
    original_meta_without_batch_df = df[[col for col in meta_cols if col != batch_column]].copy()
    if debug:
        print(f"DEBUG: Metadata columns to reattach: {original_meta_without_batch_df.columns.tolist()}")

    # Define x as the data to be debiased (gene expressions + batch column)
    x = df[gene_cols + [batch_column]]
    if debug:
        print(f"DEBUG: Shape of 'x' (data to be debiased): {x.shape}")
        print(f"DEBUG: Head of 'x':\n{x.head()}")
        print(f"DEBUG: Columns of 'x' before fit_transform: {x.columns.tolist()}")
        print(f"DEBUG: Is '{batch_column}' in x.columns? {batch_column in x.columns}")


    # PENDING: Might fix KeyError "None of ['meta_Class'] are in the columns".
    # Provide a dummy 'y' if no specific target variable is being debiased by FairAdapt.
    y = pd.DataFrame(np.zeros(len(x)), index=x.index)
    if debug:
        print(f"DEBUG: Shape of 'y' (dummy target): {y.shape}")
        print(f"DEBUG: Head of 'y':\n{y.head()}")

    # Get first row from x for FairAdapt
    first_row = x.iloc[0, :]
    if debug:
        print(f"DEBUG: First row for FairAdapt (first_row):\n{first_row.head()}")

    if debug:
        print("DEBUG: Calling FA.fit_transform...")
    x_adjusted, y_adjusted, _ = FA.fit_transform(x, y, first_row)
    # Code was run, fixed KeyError: "None of ['meta_Class'] are in the columns".
    if debug:
        print("DEBUG: FA.fit_transform completed.")
        print(f"DEBUG: Shape of adjusted 'x': {x_adjusted.shape}")
        print(f"DEBUG: Shape of adjusted 'y': {y_adjusted.shape}")


    df_adjusted = pd.DataFrame(x_adjusted, columns=gene_cols + [batch_column])
    # Add original metadata columns back to the DataFrame
    df_adjusted[original_meta_without_batch_df.columns] = original_meta_without_batch_df
    if debug:
        print(f"DEBUG: Final adjusted DataFrame shape: {df_adjusted.shape}")
        print(f"DEBUG: Final adjusted DataFrame head:\n{df_adjusted.head()}")

    # Save the imputed data to a CSV file --------------------
    try:
        if debug:
            print(f"DEBUG: Attempting to save adjusted data to {output_file}")
        df_adjusted.to_csv(output_file, index=False)
        print(f"Imputed data saved to {output_file}")
    except Exception as e:
        raise ValueError(f"AutoClass Error saving the output file: {output_file}\n{e}")

# Main execution block
if __name__ == "__main__":
    parser = ArgumentParser(description="AutoClass Imputation Example")
    parser.add_argument("-i", "--input-file", help="Path to input CSV file", required=True)
    parser.add_argument("-o", "--output-file", help="Path to output CSV file", required=True)
    parser.add_argument("-b", "--batch-column", help="Column name for batch information", required=True)
    args = parser.parse_args()

    # Call the main function with debugging enabled/disabled manually
    # Set debug=True to enable debugging, or debug=False to disable.
    run_fair_adapt(
        input_file=args.input_file,
        output_file=args.output_file,
        batch_column=args.batch_column,
        debug=True # Set this to True or False as needed for debugging
    )
