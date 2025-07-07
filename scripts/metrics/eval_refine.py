from prepdata/read_from_h5 import read_expression_data

def print_basic_stats(df):
    """
    Prints basic statistics of the DataFrame.
    """
    print("Basic Statistics:")
    print(df.describe())
    print("\nMissing Values:")
    print(df.isnull().sum())
    print("\nData Types:")
    print(df.dtypes)


def make_dataframe_from_h5(h5_path, dataset_name):
    """
    Reads expression data from an HDF5 file and returns it as a pandas DataFrame.
    
    Parameters:
    - h5_path: Path to the HDF5 file.
    - dataset_name: Name of the dataset within the HDF5 file.
    
    Returns:
    - DataFrame containing the expression data.
    """
    df = read_expression_data(h5_path, dataset_name)
    return df

if __name__ == "__main__":
    # Example usage
    h5_path = '/data/refinebio.h5'
    dataset_name = 'gse20194'
    df = make_dataframe_from_h5(h5_path, dataset_name)
    print_basic_stats(df)
