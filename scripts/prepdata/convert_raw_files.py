# /scripts/prepdata/convert_raw_files.py

import os
import pandas as pd
import gzip
import warnings
import sys
from pathlib import Path
import re
import argparse
import time
import shutil
import csv
import zipfile

# Suppress pandas warnings for cleaner output
warnings.filterwarnings('ignore', category=pd.errors.DtypeWarning)


def smart_read_dataframe(file_path: Path, debug: bool = False, **kwargs) -> tuple:
    """
    Attempt to read a dataframe using pandas' inference. If it fails,
    run diagnostics to suggest the cause of the error.
    
    Returns a (dataframe, status_string) tuple. On failure, the dataframe
    is None and the string contains diagnostic information.
    """
    # First, try common delimiters explicitly for TSV/CSV files
    for delimiter in ['\t', ',', ';']:
        try:
            df = pd.read_csv(
                file_path,
                delimiter=delimiter,
                compression='infer',
                **kwargs
            )
            if df is not None and not df.empty and len(df.columns) > 1:
                return df, f"Success (delimiter: '{delimiter}')"
        except Exception as e:
            if debug:
                print(f"DEBUG: Failed with delimiter '{delimiter}'. Error: {e}")
            continue
    
    # Fallback: try letting pandas infer everything.
    try:
        df = pd.read_csv(
            file_path,
            delimiter=None,
            engine='python',
            compression='infer',
            **kwargs
        )
        if df is not None and not df.empty:
            return df, "Success (Pandas auto-inference)"
    except Exception as e:
        if debug:
            print(f"DEBUG: Initial pandas read failed. Error: {e}")
        pass

    # If the first attempt fails, run diagnostics.
    diagnostics = []
    
    # Check 1: Validate file existence and size.
    if not file_path.exists():
        return None, "Failed: File does not exist."
    if file_path.stat().st_size == 0:
        return None, "Failed: File is empty."

    # Check 2: Check for compression mismatch (magic bytes vs. extension).
    is_gzip, is_zip = False, False
    try:
        with open(file_path, 'rb') as f:
            header = f.read(4)
            is_gzip = header.startswith(b'\x1f\x8b')
            is_zip = header.startswith(b'\x50\x4b\x03\x04')
            
        ext = file_path.suffix.lower()
        if ext == '.gz' and not is_gzip:
            diagnostics.append("File has '.gz' extension but is not a valid gzip file.")
        elif ext != '.gz' and is_gzip:
            diagnostics.append("File appears to be gzipped but lacks '.gz' extension.")
        
        if ext == '.zip' and not is_zip:
            diagnostics.append("File has '.zip' extension but is not a valid zip archive.")
        elif ext != '.zip' and is_zip:
            diagnostics.append("File appears to be a zip archive but lacks '.zip' extension.")

    except IOError as e:
        return None, f"Failed: Could not read file for diagnostics. Error: {e}"

    # Check 3: Sniff the first line for a potential delimiter.
    try:
        open_func = gzip.open if is_gzip else open
        
        if is_zip:
            diagnostics.append("Delimiter check skipped for zip archives.")
        else:
            with open_func(file_path, 'rt', encoding='utf-8', errors='ignore') as f:
                first_line = f.readline()
                sniffer = csv.Sniffer()
                dialect = sniffer.sniff(first_line, delimiters=',\t;| ')
                diagnostics.append(f"Detected potential delimiter: '{dialect.delimiter}'.")
    except (csv.Error, StopIteration):
        diagnostics.append("Could not detect a common delimiter in the first line.")
    except Exception as e:
        if debug:
            print(f"DEBUG: Delimiter check failed. Error: {e}")
        diagnostics.append("Could not read first line (possible encoding or corruption issue).")

    # Format the final diagnostic report.
    if not diagnostics:
        diagnostics.append("Pandas parser failed for an unknown reason. The file may be malformed.")

    return None, f"Failed. Diagnostics: {'; '.join(diagnostics)}"
    
# === DATA PROCESSING LOGIC ===

def _find_dataset_file(folder_path: Path, dataset_id: str, prefix: str) -> Path | None:
    """Find a dataset file with a given prefix, trying common extensions and patterns."""
    for ext in ['.csv', '.tsv']:
        # Try exact match first: e.g., expression_GSE12345.csv
        exact_path = folder_path / f"{prefix}{dataset_id.upper()}{ext}"
        if exact_path.exists():
            return exact_path
        # Fallback to pattern matching: e.g., expression_data_GSE12345_raw.tsv
        for candidate in folder_path.glob(f"{prefix}*{ext}"):
            if dataset_id.upper().replace('_', '') in candidate.stem.upper().replace('_', ''):
                return candidate
    return None

def find_dataset_files(raw_folder_path: Path, dataset_id: str) -> tuple[Path | None, Path | None]:
    """Find expression and metadata files for a given dataset ID."""
    expression_file = _find_dataset_file(raw_folder_path, dataset_id, "expression_")
    meta_file = _find_dataset_file(raw_folder_path, dataset_id, "meta_")
    return expression_file, meta_file

def _should_transpose(expr_df: pd.DataFrame, meta_df: pd.DataFrame | None) -> tuple[bool, str]:
    """Apply heuristics to determine if the expression dataframe needs transposing."""
    rows, cols = expr_df.shape
    
    # Heuristic 1: Compare dimensions to sample count from metadata
    if meta_df is not None:
        expected_samples = len(meta_df)
        # If columns are a much closer match to sample count than rows, transpose.
        if abs(cols - expected_samples) < abs(rows - expected_samples) * 0.5:
            return True, f"columns ({cols}) closer to expected samples ({expected_samples}) than rows ({rows})"
    
    # Heuristic 2: Genomics convention (genes >> samples)
    if rows > cols * 10:
        return True, f"many more rows ({rows}) than columns ({cols}), likely genes in rows"
    if cols > rows * 10:
        return False, f"many more columns ({cols}) than rows ({rows}), likely genes in columns"

    # Heuristic 3: Detect sample-like IDs in column names
    sample_like_cols = [c for c in expr_df.columns if isinstance(c, str) and c.startswith(('GSM', 'TCGA', 'Sample'))]
    if len(sample_like_cols) / cols > 0.5: # If >50% of columns look like samples
        return True, f"detected {len(sample_like_cols)} sample-like column names"
        
    # Heuristic 4: Detect gene-like IDs in the first column's content
    first_col_content = expr_df.iloc[:10, 0].astype(str).str.upper()
    if first_col_content.str.contains('ENSG|GENE|PROBE').sum() > 3:
        return True, f"first column '{expr_df.columns[0]}' contains gene-like identifiers"
        
    # Default Fallback: Assume genes are in rows if there are more rows than columns
    if rows > cols:
        return True, f"fallback: more rows ({rows}) than columns ({cols})"
        
    return False, "fallback: not transposing"

def process_dataset(raw_folder_path: Path, dataset_id: str, output_base_dir: Path, debug: bool):
    """Process a single dataset: combine expression and metadata, save to the output directory."""
    print(f"\n🔄 Processing dataset: {dataset_id}")
    
    try:
        # 1. Find and read files (this is the first validation step)
        expression_file, meta_file = find_dataset_files(raw_folder_path, dataset_id)
        if not expression_file or not meta_file:
            print("   ❌ Validation failed: Expression or metadata file not found. Skipping.")
            return False
            
        expr_df, _ = smart_read_dataframe(expression_file)
        meta_df, _ = smart_read_dataframe(meta_file)
        if expr_df is None or meta_df is None:
            print("   ❌ Validation failed: Could not read expression or metadata file. Skipping.")
            return False

        # 2. Transpose expression data if necessary
        needs_transpose, reason = _should_transpose(expr_df, meta_df)
        print(f"   🔍 Transpose decision: {'YES' if needs_transpose else 'NO'} - {reason}")
        if needs_transpose:
            # Identify gene info columns vs sample data columns
            gene_info_cols = []
            for col in expr_df.columns:
                col_lower = str(col).lower()
                if any(keyword in col_lower for keyword in ['dataset', 'entrez', 'gene', 'symbol', 'ensembl', 'chromosome', 'biotype']):
                    gene_info_cols.append(col)
                else:
                    break  # Stop at first non-gene-info column
            
            if gene_info_cols:
                # Find the best gene identifier column to use as column names after transpose
                gene_id_col = None
                for preferred in ['HGNC_Symbol', 'Entrez_Gene_ID', 'Ensembl_Gene_ID']:
                    if preferred in gene_info_cols:
                        gene_id_col = preferred
                        break
                if gene_id_col is None:
                    gene_id_col = gene_info_cols[1] if len(gene_info_cols) > 1 else gene_info_cols[0]
                
                # Create gene identifiers for column names
                gene_ids = expr_df[gene_id_col].astype(str)
                
                # Separate expression data and transpose
                expr_data = expr_df.drop(columns=gene_info_cols)
                expr_df = expr_data.T
                expr_df.columns = gene_ids
                expr_df.index.name = 'Sample_ID'
            else:
                # Fallback: use first column as gene index
                expr_df = expr_df.set_index(expr_df.columns[0]).T
        else:
            # Assume first column is sample ID
            expr_df = expr_df.set_index(expr_df.columns[0])
        
        # 3. Prepare metadata
        # Find and set sample ID column as index - prioritize 'Sample_ID' specifically
        sample_id_col = None
        for col in meta_df.columns:
            col_lower = str(col).lower()
            if col_lower == 'sample_id':
                sample_id_col = col
                break
            elif 'sample' in col_lower and 'id' in col_lower:
                sample_id_col = col
                break
        
        # If no sample_id column found, look for columns with GSM-like values
        if sample_id_col is None:
            for col in meta_df.columns:
                if meta_df[col].astype(str).str.contains('GSM|TCGA|Sample').any():
                    sample_id_col = col
                    break
        
        if sample_id_col is None:
            # Last resort: use second column if first is dataset_id
            if len(meta_df.columns) > 1 and 'dataset' in str(meta_df.columns[0]).lower():
                sample_id_col = meta_df.columns[1]
            else:
                sample_id_col = meta_df.columns[0]
            
        if debug:
            print(f"DEBUG: Using metadata column '{sample_id_col}' as sample ID")
            
        meta_df = meta_df.set_index(sample_id_col)
        # Prefix columns with 'meta_'
        meta_df = meta_df.add_prefix('meta_')

        lower_to_full_case_columns = {col.lower():col for col in meta_df.columns}

        # Standardize er status column name
        other_er_columns = ['meta_er_ihc', 'meta_er', 'meta_er_status', 'meta_er_status_by_ihc', 'meta_er_status_ihc', 'meta_esr1_status']
        for er_column in other_er_columns:
            if er_column in lower_to_full_case_columns:
                er_column = lower_to_full_case_columns[er_column]
                print(f"Using {er_column} as ER status column")
                meta_df = meta_df.rename(columns={er_column: 'meta_er_status'})
                break

        if "meta_er_status" not in meta_df.columns.tolist():
            print(f"   ❌ ER status column not found. ")
            print(f"DEBUG: Available columns: {meta_df.columns.tolist()}")
        else:
            print(f"Unique er values: ")
            print(meta_df['meta_er_status'].value_counts())

        # Standardize pr status column name
        other_pr_columns = ['meta_pr', 'meta_pr_status', 'meta_pr_status_ihc', 'meta_pr_ihc', 'meta_prihc', 'meta_per_status_by_ihc']
        for pr_column in other_pr_columns:
            if pr_column in lower_to_full_case_columns:
                pr_column = lower_to_full_case_columns[pr_column]
                print(f"Using {pr_column} as PR status column")
                meta_df = meta_df.rename(columns={pr_column: 'meta_pr_status'})
                break

        if "meta_pr_status" not in meta_df.columns.tolist():
            print(f"   ❌ PR status column not found. ")
            print(f"DEBUG: Available columns: {meta_df.columns.tolist()}")
        else:
            print(f"Unique pr values: ")
            print(meta_df['meta_pr_status'].value_counts())

        # Standardize her2 status column name
        other_her2_columns = ['meta_her2', 'meta_her_2', 'meta_her2_status', 'meta_her_2_status', 'meta_her2_status_by_ihc']
        for her2_column in other_her2_columns:
            if her2_column in lower_to_full_case_columns:
                her2_column = lower_to_full_case_columns[her2_column]
                print(f"Using {her2_column} as HER2 status column")
                meta_df = meta_df.rename(columns={her2_column: 'meta_her_2_status'})
                break

        if "meta_her2_status" not in meta_df.columns.tolist():
            print(f"   ❌ HER2 status column not found. ")
            print(f"DEBUG: Available columns: {meta_df.columns.tolist()}")
        else:
            print(f"Unique her2 values: ")
            print(meta_df['meta_her2_status'].value_counts())

        # Map column values to 0 for negative and 1 for positive.
        expected_cols = ['meta_er_status', 'meta_pr_status', 'meta_her2_status']
        cols_to_convert = [col for col in expected_cols if col in meta_df.columns]

        # Lowercase all values in relevant columns
        meta_df[cols_to_convert] = meta_df[cols_to_convert].apply(lambda col: col.str.lower())

        def status_to_binary(val):
            if pd.isnull(val):
                return np.nan 
            # Case 1: Already numeric
            if isinstance(val, (int, float)):
                if val == 0.0 or val == 0:
                    return 0
                if val == 1.0 or val == 1:
                    return 1

            # Case 2: It's a string
            if isinstanace(val, str):

                val = val.strip().lower()

                positive_vals = {'positive', 'P', 'pos', 'pos-low', '1', 'er+', 'he+', 'pr+'}
                negative_vals = {'negative', 'N', 'neg', '0', 'er-', 'he-', 'pr-'}

                for pos in positive_vals:
                    if pos in val:
                        return 1
                for neg in negative_vals:
                    if neg in val:
                        return 0

            return None

        df[cols_to_convert] = df[cols_to_convert].applymap(status_to_binary)

        print("Number of values that couldn't be classified: ", df[cols_to_convert].isnull().sum())

        # 4. Align and combine (the final and most critical validation)
        common_samples = expr_df.index.intersection(meta_df.index)
        if len(common_samples) == 0:
            print(f"   ❌ Validation failed: No common samples found between files. Skipping.")
            if debug:
                # DEBUG: Show sample ID examples to help diagnose mismatches.
                expr_samples = sorted(list(expr_df.index))[:5]
                meta_samples = sorted(list(meta_df.index))[:5]
                print(f"DEBUG: Expression sample IDs: {expr_samples}")
                print(f"DEBUG: Metadata sample IDs: {meta_samples}")
                print(f"DEBUG: Expression index type: {type(expr_df.index[0]) if len(expr_df.index) > 0 else 'empty'}")
                print(f"DEBUG: Metadata index type: {type(meta_df.index[0]) if len(meta_df.index) > 0 else 'empty'}")
            return False
            
        print(f"   🔗 Found {len(common_samples)} common samples.")
        combined_df = pd.concat([expr_df.loc[common_samples], meta_df.loc[common_samples]], axis=1)

        # 5. Save result
        output_dir = output_base_dir / dataset_id.lower()
        output_dir.mkdir(parents=True, exist_ok=True)
        output_file = output_dir / "unadjusted.csv"
        combined_df.to_csv(output_file)
        
        print(f"   ✅ Successfully saved: {output_file} ({output_file.stat().st_size:,} bytes)")
        return True

    except Exception as e:
        print(f"   ❌ Error processing {dataset_id}: {e}")
        if debug:
            import traceback
            traceback.print_exc()
        return False

def scan_for_datasets(raw_data_dir: Path) -> list[dict]:
    """Scan the raw data directory for valid dataset folders."""
    datasets = []
    if not raw_data_dir.exists():
        print(f"Warning: Raw data directory not found at {raw_data_dir}")
        return datasets
    for item_path in raw_data_dir.iterdir():
        if item_path.is_dir():
            expression_file, meta_file = find_dataset_files(item_path, item_path.name)
            if expression_file and meta_file:
                datasets.append({'folder': item_path, 'dataset_id': item_path.name})
    return datasets

# --- MAIN EXECUTION BLOCK ---
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process raw genomics data into combined datasets.")
    # Add arguments for directories and debug flag
    parser.add_argument('--raw-dir', type=Path, required=True,
                        help='Directory containing the raw dataset folders.')
    parser.add_argument('--target-dir', type=Path, required=True,
                        help='Directory to save the processed output files.')
    parser.add_argument('--debug', action='store_true', help='Enable detailed debug output.')
    args = parser.parse_args()
    
    print("="*80)
    print(f"PROCESSING RAW DATASETS -> {args.target_dir}", flush=True)
    print("="*80)

    datasets_info = scan_for_datasets(args.raw_dir)
    print(f"\nFound {len(datasets_info)} datasets to process in {args.raw_dir}:")
    for info in datasets_info:
        print(f"  - {info['dataset_id']}")
    
    # Process each dataset
    successful_count = 0
    for info in datasets_info:
        if process_dataset(info['folder'], info['dataset_id'], args.target_dir, args.debug):
            successful_count += 1

    print("\n" + "="*80)
    print("PROCESSING SUMMARY")
    print("="*80)
    print(f"✅ Successfully processed: {successful_count}/{len(datasets_info)} datasets")
    print(f"❌ Failed to process: {len(datasets_info) - successful_count}/{len(datasets_info)} datasets")
    print(f"📁 Output directory: {args.target_dir}")
