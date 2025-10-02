# combine_datasets.py - Simplified version with gene ID conversion moved upstream

import argparse
import os
import pandas as pd

def print_now(*args, **kwargs):
    """Prints a message to the console with flushing to ensure immediate output."""
    print(*args, flush=True, **kwargs)

def get_common_genes(df1, df2):
    """Finds the intersection of gene columns between two dataframes."""
    genes1 = [col for col in df1.columns if not col.startswith('meta_') and col != "Sample_ID"]
    genes2 = [col for col in df2.columns if not col.startswith('meta_') and col != "Sample_ID"]
    common_genes = list(set(genes1) & set(genes2))
    return common_genes, genes1, genes2

def create_result(success=True, error=None, **kwargs):
    """Create a standardized result dictionary."""
    base_result = {
        'success': success,
        'error': error,
        'dataset1_samples': 0,
        'dataset2_samples': 0,
        'combined_samples': 0,
        'common_genes': 0,
        'genes1_count': 0,
        'genes2_count': 0
    }
    base_result.update(kwargs)
    return base_result

def combine_datasets(input1_path, input2_path, output_path):
    """Combine two datasets and return structured results."""
    try:
        # Load datasets
        df1 = pd.read_csv(input1_path, low_memory=False)
        df2 = pd.read_csv(input2_path, low_memory=False)
        
        # Validate required columns
        missing_er_status = [i for i, df in enumerate([df1, df2], 1) if 'meta_er_status' not in df.columns]
        if missing_er_status:
            return create_result(
                success=False, 
                error='meta_er_status_missing',
                error_details={'missing_datasets': missing_er_status},
                dataset1_samples=df1.shape[0],
                dataset2_samples=df2.shape[0]
            )
        
        # Get common genes
        common_genes, genes1, genes2 = get_common_genes(df1, df2)
        if not common_genes:
            return create_result(
                success=False,
                error='no_common_genes', 
                dataset1_samples=df1.shape[0],
                dataset2_samples=df2.shape[0],
                genes1_count=len(genes1),
                genes2_count=len(genes2)
            )
        
        # Add source identifiers and combine
        df1['meta_source'] = os.path.basename(os.path.dirname(input1_path))
        df2['meta_source'] = os.path.basename(os.path.dirname(input2_path))
        
        keep_cols = common_genes + ['meta_er_status', 'meta_source']
        df_combined = pd.concat([df1[keep_cols], df2[keep_cols]], ignore_index=True)
        
        # Save result
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        df_combined.to_csv(output_path, index=False)
        
        return create_result(
            dataset1_samples=df1.shape[0],
            dataset2_samples=df2.shape[0],
            combined_samples=df_combined.shape[0],
            common_genes=len(common_genes),
            genes1_count=len(genes1),
            genes2_count=len(genes2),
            output_path=output_path
        )
        
    except Exception as e:
        return create_result(
            success=False,
            error='exception',
            error_details={'exception_type': type(e).__name__, 'exception_message': str(e)}
        )

def format_result(result):
    """Format combination result for display."""
    if result['success']:
        msg = f"SUCCESS: {result['combined_samples']} samples, {result['common_genes']} genes"
        
        warnings = []
        if result['common_genes'] < 10:
            warnings.append(f"very few genes ({result['common_genes']})")
        elif result['common_genes'] < 100:
            warnings.append(f"low gene count ({result['common_genes']})")
        if result['combined_samples'] < 20:
            warnings.append(f"small sample size ({result['combined_samples']})")
        
        return msg + (f" ⚠️ {', '.join(warnings)}" if warnings else "")
    
    # Format error messages
    error_type = result['error']
    details = result.get('error_details', {})
    
    error_messages = {
        'meta_er_status_missing': lambda d: f"'meta_er_status' column missing in dataset(s) {', '.join(map(str, d.get('missing_datasets', [])))}",
        'no_common_genes': lambda d: "No common genes found between datasets",
        'exception': lambda d: f"{d.get('exception_type', 'Exception')}: {d.get('exception_message', 'Unknown error')}"
    }
    
    error_msg = error_messages.get(error_type, lambda d: str(error_type))(details)
    return f"ERROR: {error_msg}"

def main():
    """Main function to parse arguments and orchestrate the dataset combination."""
    parser = argparse.ArgumentParser(description="Combine two gene expression datasets based on common genes.")
    parser.add_argument('--input1', required=True, help='Path to the first data file.')
    parser.add_argument('--input2', required=True, help='Path to the second data file.')
    parser.add_argument('--output', required=True, help='Path for the combined output CSV file.')
    parser.add_argument('--debug', action='store_true', help='Enable detailed debug output.')

    args = parser.parse_args()

    result = combine_datasets(args.input1, args.input2, args.output)
    print_now(format_result(result))
    
    return 0 if result['success'] else 1

if __name__ == "__main__":
    main()