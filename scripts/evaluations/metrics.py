from sklearn.metrics import make_scorer, roc_auc_score, accuracy_score, log_loss, mutual_info_score
import numpy as np
from sklearn.feature_selection import mutual_info_classif
from sklearn.model_selection import cross_val_score

### METRICS ###

def mutual_info_shannons(y_true, y_pred):
    """Calculate mutual information in shannons (bits)."""
    return mutual_info_score(y_true, y_pred) / np.log(2)

def mutual_info_proba_classif(y_true, y_pred_proba):
    """Calculate MI for categorical true labels and continuous predictions."""
    mi_value = mutual_info_classif(y_pred_proba.reshape(-1, 1), y_true, random_state=42)[0]
    return mi_value / np.log(2)


def mutual_info_proba_sum(y_true, y_pred_proba, classes=[]):
    """
    Calculate mutual information by summing the probabilities across samples, for each true class.
    This results in a contingency table (2x2 for binary classification), which is then used to compute mutual information.
    This is similar to mutual_info_score but sums p(y_pred_proba | y_true) across samples, instead of 1-0 encoding the probabilities.
    """
    # For binary classification, y_pred_proba is a 1D array: [n_samples]
    # For multiclass classification, y_pred_proba is a 2D array: [n_samples, n_classes]
    binary = len(y_pred_proba.shape) == 1

    if not classes:
        classes = sorted(np.unique(y_true).tolist())

    if len(y_pred_proba.shape) == 1: # Binary classification
        y_pred_proba = np.column_stack((y_pred_proba, 1 - y_pred_proba))

    contingency_table = {cls: np.zeros(len(classes)) for cls in classes}
    for true_class, proba in zip(y_true, y_pred_proba):
        contingency_table[true_class] += proba

    # Scale to avoid precision loss when converting to int for contingency score
    contingency_array = np.array(list(contingency_table.values())) * 10**6
    return mutual_info_score(None, None, contingency=contingency_array) / np.log(2)

def one_minus_log_loss(y_true, y_pred_proba):
    """
    Calculate one minus log loss.
    This is useful for maximizing the score, as log loss is minimized.
    """
    return 1 - log_loss(y_true, y_pred_proba)


def determinant_based_mutual_information(y_true, y_pred_proba):
    """
    Calculate mutual information using det(P), where P is the joint probability distribution of y_true and y_pred_proba.
    """
    classes = sorted(np.unique(y_true).tolist())

    # Convert y_true to one-hot encoding
    one_hot = np.zeros((len(y_true), len(classes)))
    for i, cls in enumerate(classes):
        one_hot[:, i] = (y_true == cls).astype(int)
    y_true = one_hot

    y_pred_proba = np.asarray(y_pred_proba)
    
    binary = len(y_pred_proba.shape) == 1
    if binary:
        # For binary classification, y_pred_proba is a 1D array: [n_samples]
        # For multiclass classification, y_pred_proba is a 2D array: [n_samples, n_classes]
        y_pred_proba = np.column_stack((y_pred_proba, 1 - y_pred_proba))

    contingency_table = y_pred_proba.T @ y_true
    rel_freq_table = contingency_table / contingency_table.sum()

    # Calculate the determinant of the contingency table
    det = np.linalg.det(rel_freq_table + 1e-10)  # Add a small value to avoid singular matrix issues
    return abs(det) * 4 # Max is 0.25 (max of x*1-x) because the probabilities must sum to 1. [0.5, 0] x [0, 0.5] = 0.25

