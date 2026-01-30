#!/usr/bin/env python3

import argparse
import pandas as pd
import re
import numpy as np

def load_data(input_csv):
    return pd.read_csv(input_csv)


def load_metadata(metadata_file):
    return pd.read_csv(metadata_file)


def standardize_columns(df, metadata):
    """
    Standardize columns using final_name, merging dataset-specific columns.
    If no source columns exist in df, create final_name column as NaN.
    """
    from collections import defaultdict

    dataset_cols = [c for c in metadata.columns if c != "final_name"]

    for _, row in metadata.iterrows():
        final_name = row["final_name"]
        # List all dataset-specific columns for this final_name
        sources = [row[c] for c in dataset_cols if pd.notna(row[c])]

        # Keep only columns that exist in the df
        existing_sources = [c for c in sources if c in df.columns]

        if existing_sources:
            if len(existing_sources) == 1:
                # Only one exists → rename to final_name
                src = existing_sources[0]
                if src != final_name:
                    df.rename(columns={src: final_name}, inplace=True)
            else:
                # Multiple exist → merge with bfill
                df[final_name] = df[existing_sources].bfill(axis=1).iloc[:, 0]
                df.drop(columns=existing_sources, inplace=True)
                print(f"🔀 Merged {existing_sources} → {final_name}")
        else:
            # None exist → create empty column
            df[final_name] = np.nan
            print(f"⚠️ Created empty column {final_name} (no sources present)")

    return df

def map_status_to_binary(df, columns):
    """
    Convert metadata columns to 0/1.
    """
    # Keep only existing columns
    cols_to_convert = [col for col in columns if col in df.columns]
    if not cols_to_convert:
        return df

    # Lowercase all string values
    df[cols_to_convert] = df[cols_to_convert].apply(lambda col: col.astype(str).str.lower())

    def status_to_binary(val):
        if pd.isnull(val):
            return np.nan

        val = str(val).strip().lower()

        # Try numeric first
        try:
            num_val = float(val)
            if num_val == 0:
                return 0
            elif num_val in [1, 2, 3]:
                return 1
        except ValueError:
            pass

        positive_vals = {'male','positive', 'yes'}
        negative_vals = {'female','negative', 'no'}

        for pos in positive_vals:
            if pos in val:
                return 1
        for neg in negative_vals:
            if neg in val:
                return 0

        return np.nan

    for col in cols_to_convert:
        df[col] = df[col].map(status_to_binary)

        # Report unclassified values
        unclassified = df[col][df[col].isnull()]
        if not unclassified.empty:
            print(f"❓ Unclassified values in {col}:")
            print(unclassified.value_counts(dropna=False))

        # Show post-conversion unique values
        print(f"✅ Post-conversion unique values in {col}:")
        print(df[col].value_counts(dropna=False))

    return df

def map_column_with_regex(df, column_name, patterns_to_values):
    """
    Map values in a column to standardized numeric or categorical values
    using regex patterns.

    Parameters:
        df (pd.DataFrame): The dataframe
        column_name (str): The column to transform
        patterns_to_values (list of tuples): [(regex_pattern, mapped_value), ...]

    Returns:
        pd.Series: Transformed column
    """
    if column_name not in df.columns:
        print(f"⚠️ Column {column_name} not found, skipping regex mapping")
        return df

    def map_value(val):
        if pd.isnull(val):
            return np.nan
        val_str = str(val).lower().strip()
        for pattern, mapped_val in patterns_to_values:
            if re.search(pattern, val_str):
                return mapped_val
        return np.nan  # fallback if no pattern matched

    df[column_name] = df[column_name].apply(map_value)

    # Optional: report unmapped values
    unmapped = df[column_name][df[column_name].isnull()]
    if not unmapped.empty:
        print(f"❓ Unmapped values in {column_name}:")
        print(unmapped.value_counts(dropna=False))

    return df

def one_hot_encode_age(
    df,
    age_col="meta_age_at_diagnosis",
    bins=(0, 40, 50, 60, 70, 200),
    labels=("lt40", "40_49", "50_59", "60_69", "ge70"),
    drop_original=True
):
    """
    Bin age_at_diagnosis and one-hot encode the bins.
    """
    if age_col not in df.columns:
        print(f"⚠️ {age_col} not found, skipping age one-hot encoding")
        return df

    # Force numeric
    df[age_col] = pd.to_numeric(df[age_col], errors="coerce")

    # Bin ages
    age_bins = pd.cut(
        df[age_col],
        bins=bins,
        labels=labels,
        right=False
    )

    # One-hot encode
    age_dummies = pd.get_dummies(
        age_bins,
        prefix=age_col,
        dummy_na=False
    )

    df = pd.concat([df, age_dummies], axis=1)

    if drop_original:
        df = df.drop(columns=[age_col])

    print("✅ Age one-hot encoding complete:")
    print(age_dummies.sum())

    return df

def write_output(df, output_path):
    df.to_csv(output_path, index=False)
    print(f">>> Standardized file written to: {output_path}")

def normalize_colname(c):
    return c.strip().lower()

def main():
    parser = argparse.ArgumentParser(
        description="Standardize column names and map categorical status to numeric"
    )
    parser.add_argument("--input", required=True,
                        help="Input CSV file")
    parser.add_argument("--metadata_file", required=True,
                        help="CSV file with column name mappings")
    parser.add_argument("--output", required=True,
                        help="Output CSV file")
    args = parser.parse_args()

    df = load_data(args.input)
    metadata = load_metadata(args.metadata_file)

    print("After loading:")
    meta_cols = [c for c in df.columns if c.startswith("meta_")]
    print(meta_cols)

    df.columns = [normalize_colname(c) for c in df.columns]
    metadata.columns = [normalize_colname(c) for c in metadata.columns]

    for col in metadata.columns:
        if col != "final_name":
            metadata[col] = metadata[col].astype(str).str.strip().str.lower()

    # Step 1: standardize column names
    df = standardize_columns(df, metadata)

    print("After standardizing")
    meta_cols = [c for c in df.columns if c.startswith("meta_")]

    # Print just those columns
    print(df[meta_cols].head())

    # Step 2: map tricky columns with regex
    df = map_column_with_regex(
        df, 
        column_name='meta_menopause_status',
        patterns_to_values=[
            (r"pre", 0),
            (r"post", 1),
        ]
    )

    df = map_column_with_regex(
        df, 
        column_name='meta_her2_status',
        patterns_to_values=[
            (r"equivocal", 1),
            (r"her2-|negative", 0),
            (r"her2\+|positive", 2),
        ]
    )

    df = map_column_with_regex(
        df, 
        column_name='meta_histological_type',
        patterns_to_values=[
            (r"infiltrating ductal|ductal", 0),
            (r"infiltrating lobular|lobular", 1),
            (r"medullary", 2),
            (r"metaplastic", 3),
            (r"mixed", 4),
            (r"mucinous", 5),
            (r"other|nos", 6),
            (r"tubular", 7)
        ]
    )

    # Unique GSE62944_TUMOR histological_type values:
    # Infiltrating Carcinoma NOS
    # Infiltrating Ductal Carcinoma
    # Infiltrating Lobular Carcinoma
    # Medullary Carcinoma
    # Metaplastic Carcinoma
    # Mixed Histology (please specify)
    # Mucinous Carcinoma
    # NA
    # Other  specify

    # Unique METABRIC histological_subtype values:
    # Ductal/NST
    # Lobular
    # Medullary
    # Metaplastic
    # Mixed
    # Mucinous
    # NA
    # Other
    # Tubular/ cribriform

    # Step 2: one-hot encode age at diagnosis
    df = one_hot_encode_age(
        df,
        age_col="meta_age_at_diagnosis",
        bins=(0, 50, 70, 200),
        labels=("lt50", "50_69", "ge70"),
        drop_original=False 
    )


    # Step 3: convert target metadata columns to 0/1
    status_cols = ['meta_sex', 'meta_chemotherapy']
    df = map_status_to_binary(df, status_cols)

    # Step 4: write final output
    write_output(df, args.output)


if __name__ == "__main__":
    main()
