"""
Rank-In: Python implementation for integrative analysis across microarray and RNA-seq.

Reference:
    Tang K, Ji X, Zhou M, Deng Z, Huang Y, Zheng G, Cao Z.
    "Rank-in: enabling integrative analysis across microarray and RNA-seq for cancer."
    Nucleic Acids Research. 2021;49(17):e99. doi:10.1093/nar/gkab554

Algorithm (Section 2.1 of the paper):
  1. For each expression profile, rank genes from lowest to highest intensity.
  2. Partition ranked genes into N_GROUPS (100) equal-sized bins.
  3. Within each bin, weight the ranks by the linear slope of raw intensities
     across the bin (captures the "steepness" of each expression stratum).
  4. Concatenate all weighted-rank profiles into a combined matrix.
  5. Apply SVD to remove the first K (default=1) singular vectors, which
     capture platform-specific batch variation.
"""

import numpy as np


N_GROUPS = 100   # Number of rank bins (fixed per paper)
N_SVD = 1        # Number of SVD components to remove (batch PCs)


def _rank_and_weight(profile: np.ndarray) -> np.ndarray:
    """
    Convert a single expression profile to its weighted rank representation.

    Parameters
    ----------
    profile : 1-D array of shape (n_genes,)
        Raw expression values for one sample.

    Returns
    -------
    weighted : 1-D array of shape (n_genes,)
        Weighted rank values in original gene order.
    """
    n = len(profile)
    # Step 1: rank (1-based, ties broken by first occurrence)
    order = np.argsort(profile, kind="stable")       # gene indices sorted asc
    ranks = np.empty(n, dtype=float)
    ranks[order] = np.arange(1, n + 1, dtype=float)

    # Step 2: split sorted indices into N_GROUPS bins
    bin_edges = np.array_split(order, N_GROUPS)

    weighted = np.empty(n, dtype=float)
    for bin_indices in bin_edges:
        if len(bin_indices) == 0:
            continue
        bin_vals = profile[bin_indices]
        bin_ranks = ranks[bin_indices]

        # Step 3: slope of raw intensities across the sorted bin positions
        x = np.arange(len(bin_indices), dtype=float)
        if len(bin_indices) > 1:
            # Linear regression slope
            x_c = x - x.mean()
            slope = np.dot(x_c, bin_vals) / np.dot(x_c, x_c)
        else:
            slope = 1.0  # single-gene bin: neutral weight

        # Weight = rank * |slope| (slope encodes how "spread" this stratum is)
        weight = abs(slope) if slope != 0 else 1e-8
        weighted[bin_indices] = bin_ranks * weight

    return weighted


def rank_in(matrix: np.ndarray, n_svd: int = N_SVD) -> np.ndarray:
    """
    Apply the Rank-In algorithm to a gene-expression matrix.

    Parameters
    ----------
    matrix : np.ndarray of shape (n_genes, n_samples)
        Raw expression matrix. Rows = genes, columns = samples.
    n_svd : int, optional
        Number of leading singular vectors (batch PCs) to remove.
        Default = 1 per the paper's recommendation.

    Returns
    -------
    corrected : np.ndarray of shape (n_genes, n_samples)
        Batch-corrected expression matrix in the same shape as input.
    """
    n_genes, n_samples = matrix.shape

    # Step 1-3: Compute weighted-rank profile for each sample
    weighted = np.column_stack([
        _rank_and_weight(matrix[:, j]) for j in range(n_samples)
    ])  # shape: (n_genes, n_samples)

    # Step 4-5: SVD on the combined weighted-rank matrix, remove top batch PCs
    # Center rows (genes) before SVD to preserve biological differences between genes
    row_means = weighted.mean(axis=1, keepdims=True)
    centered = weighted - row_means

    U, S, Vt = np.linalg.svd(centered, full_matrices=False)

    # Reconstruct without the first n_svd components
    S_reduced = S.copy()
    S_reduced[:n_svd] = 0.0
    corrected = U @ np.diag(S_reduced) @ Vt + row_means

    return corrected


def rank_in_from_r(matrix_list, n_svd: int = N_SVD):
    """
    Entry point called from R via reticulate.

    Parameters
    ----------
    matrix_list : list of lists or numpy array
        Expression matrix passed from R (genes x samples, R convention).
    n_svd : int
        Number of SVD batch components to remove.

    Returns
    -------
    list of lists (R-compatible)
        Corrected matrix as a nested list for R to reconstruct.
    """
    mat = np.array(matrix_list, dtype=float)
    result = rank_in(mat, n_svd=int(n_svd))
    return result.tolist()
