import numpy as np
import pandas as pd
import tables
from sklearn.model_selection import StratifiedKFold, train_test_split
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.preprocessing import LabelEncoder
from sklearn.metrics import accuracy_score, roc_auc_score
import warnings

# Suppress warnings from sklearn about classes with single members in folds
warnings.filterwarnings('ignore', category=UserWarning)


def get_sample_metadata(h5_path, expression_dataset_name, metadata_dataset_name, metadata_column, debug=False):
    """
    Retrieves a metadata column and aligns it with the expression data samples.
    Handles duplicate sample IDs in the metadata by keeping the first entry.
    """
    print(f"Aligning metadata column '{metadata_column}'...")
    
    # --- Get sample order from the expression dataset ---
    with tables.open_file(h5_path, mode='r') as h5file:
        sample_names_node = h5file.get_node(f"/expression/{expression_dataset_name}", 'sample_names')
        expression_sample_names = [s.decode('utf-8') for s in sample_names_node.read()]

    # --- Read the metadata table using pandas HDFStore ---
    with pd.HDFStore(h5_path, 'r') as h5_store:
        metadata_key = f"/metadata/{metadata_dataset_name}"
        if metadata_key not in h5_store.keys():
            raise KeyError(f"Metadata key '{metadata_key}' not found in HDF5 file.")
        
        metadata_df = h5_store.get(metadata_key)
        
        if metadata_column not in metadata_df.columns:
            if debug:
                print(f"DEBUG: Available columns in metadata: {metadata_df.columns.tolist()}")
            raise KeyError(f"Column '{metadata_column}' not found in metadata.")

    # --- FIX for duplicate labels ---
    if metadata_df.index.has_duplicates:
        # Code was run, fixed "cannot reindex on an axis with duplicate labels" error. (Verified)
        metadata_df = metadata_df[~metadata_df.index.duplicated(keep='first')]

    # --- Align metadata to the expression sample order ---
    aligned_metadata = metadata_df[metadata_column].reindex(expression_sample_names)
    
    print("✓ Metadata successfully aligned.")
    return aligned_metadata


def run_platform_prediction_cv(h5_path, expression_dataset_name, y_target, n_folds=3):
    """
    Runs cross-validation to predict a target variable from expression data.

    Args:
        h5_path (str): Path to the HDF5 file.
        expression_dataset_name (str): Stem name of the expression dataset.
        y_target (pd.Series): A pandas Series with sample IDs as the index and
                              the target labels as values.
        n_folds (int): The number of folds for cross-validation.
    """
    print("\n--- Starting Platform Prediction using Cross-Validation ---")
    
    # --- 1. Preprocess Target Variable ---
    # Remove samples with missing labels
    y_target = y_target.dropna()
    if y_target.empty:
        print("🔴 ERROR: No valid labels found after dropping NaNs. Cannot proceed.")
        return

    # Filter out classes with too few samples before proceeding
    min_samples_per_class = 20
    value_counts = y_target.value_counts()
    classes_to_keep = value_counts[value_counts >= min_samples_per_class].index
    
    original_sample_count = len(y_target)
    y_target = y_target[y_target.isin(classes_to_keep)]
    
    print(f"Filtering out platforms with less than {min_samples_per_class} samples...")
    print(f"Removed {original_sample_count - len(y_target)} samples. {len(y_target)} samples remaining for analysis.")

    # --- NEW SAMPLING STRATEGY: Filter to top 10 platforms and balance classes ---
    print("Filtering to the 10 most frequent platforms...")
    top_10_platforms = y_target.value_counts().nlargest(10).index
    y_target = y_target[y_target.isin(top_10_platforms)]
    print(f"Filtered to {len(y_target)} samples from the top 10 platforms.")

    print("Balancing classes by downsampling to the smallest class size...")
    min_class_size = y_target.value_counts().min()
    print(f"Smallest class size is {min_class_size}. Downsampling all classes to this size.")

    # Use groupby and sample to balance the dataset
    y_target = (
        y_target.groupby(y_target, group_keys=False)
                .apply(lambda x: x.sample(n=min_class_size, random_state=42))
    )
    print(f"New balanced dataset size: {len(y_target)} samples ({len(y_target.unique())} classes * {min_class_size} samples).")


    # Encode string labels into integers
    le = LabelEncoder()
    y_encoded = le.fit_transform(y_target)
    print(f"Found {len(le.classes_)} unique platforms to predict in the dataset.")
    
    # Define all possible class labels that the model knows about.
    all_encoded_labels = le.transform(le.classes_)

    # --- 2. Load Expression Data (Features) ---
    print("Loading transposed expression data (samples x genes)...")
    with tables.open_file(h5_path, mode='r') as h5file:
        expression_node = h5file.get_node(f"/expression/{expression_dataset_name}", 'expression_data_transposed')
        
        # Create a mapping from sample name to its integer index for fast lookups
        all_samples_list = [s.decode('utf-8') for s in h5file.get_node(f"/expression/{expression_dataset_name}", 'sample_names').read()]
        sample_to_idx_map = {name: i for i, name in enumerate(all_samples_list)}
        
        # Use the map to get indices efficiently instead of a slow list.index() loop
        sample_indices_to_keep = [sample_to_idx_map[s] for s in y_target.index]
        
        # Use direct indexing on the PyTables node for memory-efficient row selection.
        X = expression_node[sample_indices_to_keep, :]
    print(f"Feature matrix loaded with shape: {X.shape}")

    # --- 3. Run Cross-Validation ---
    skf = StratifiedKFold(n_splits=n_folds, shuffle=True, random_state=42)
    
    accuracies = []
    auroc_scores = []
    
    for i, (train_index, test_index) in enumerate(skf.split(X, y_encoded)):
        fold_num = i + 1
        print(f"\n--- Processing Fold {fold_num}/{n_folds} ---")
        
        X_train, X_test = X[train_index], X[test_index]
        y_train, y_test = y_encoded[train_index], y_encoded[test_index]
        
        print(f"Train set size: {len(X_train)}, Test set size: {len(X_test)}")
        
        print("Training HistGradientBoostingClassifier...")
        model = HistGradientBoostingClassifier(random_state=42)
        model.fit(X_train, y_train)
        
        print("Evaluating model...")
        y_pred_proba = model.predict_proba(X_test)
        y_pred = model.predict(X_test)
        
        acc = accuracy_score(y_test, y_pred)
        accuracies.append(acc)
        
        # This handles cases where a fold's training data doesn't contain all possible classes.
        # We create a full probability matrix and insert the model's predictions
        # into the correct columns, ensuring dimensions always match.
        full_proba = np.zeros((len(y_test), len(all_encoded_labels)))
        # Code was run, fixed ValueError: 'y_true' contains labels not in parameter 'labels' (Verified)
        full_proba[:, model.classes_] = y_pred_proba
        
        auroc = roc_auc_score(y_test, full_proba, multi_class='ovr', labels=all_encoded_labels)
        auroc_scores.append(auroc)
        
        print(f"Fold {fold_num} Accuracy: {acc:.4f}")
        print(f"Fold {fold_num} AUROC:    {auroc:.4f}")

    # --- 4. Report Final Results ---
    print("\n--- Cross-Validation Results ---")
    print(f"Average Accuracy: {np.mean(accuracies):.4f} ± {np.std(accuracies):.4f}")
    print(f"Average AUROC:    {np.mean(auroc_scores):.4f} ± {np.std(auroc_scores):.4f}")


if __name__ == "__main__":
    # --- Define your dataset names ---
    H5_FILE_PATH = '/data/gold/refinebio.h5'
    EXPRESSION_NAME = 'HOMO_SAPIENS' 
    METADATA_NAME = 'metadata_HOMO_SAPIENS_indexed'
    METADATA_COLUMN_OF_INTEREST = 'refinebio_platform' 

    # 1. Get the aligned metadata first
    platforms = get_sample_metadata(
        H5_FILE_PATH,
        EXPRESSION_NAME,
        METADATA_NAME,
        METADATA_COLUMN_OF_INTEREST
    )

    print("\n--- Aligned Sample Metadata ---")
    print(platforms.head())
    print(f"\nTotal samples with metadata: {len(platforms)}")
    print(f"Samples with missing platform info: {platforms.isna().sum()}")
    print(f"Unique platforms found: {platforms.nunique()}")

    # 2. Run the machine learning cross-validation
    run_platform_prediction_cv(
        H5_FILE_PATH,
        EXPRESSION_NAME,
        y_target=platforms,
        n_folds=3
    )
