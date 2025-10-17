"""
gmm_adjust_python.py

2-component Gaussian Mixture Model (1D) with conjugate inverse-gamma variance priors
coupled across components (prior mean for each variance = other component's variance).
Uses POSTERIOR MEAN estimators for both weights (Dirichlet posterior mean) and variances
(Inverse-Gamma posterior mean). No prior on means.

Provides helpers:
- inverse_cdf_gmm
- simple_fallback
- get_gene_gmm_transform
- bimodal_normalize (with optional parallel processing)
- gmm_adjust (batch-level adjustment)
"""

import numpy as np
import pandas as pd
from scipy.stats import norm
from scipy.optimize import brentq
from multiprocessing import Pool, cpu_count
import warnings

# -----------------------
# GaussianMixture1D class
# -----------------------
# --- GMM Implementation ---
class GaussianMixture1D:
    """
    2-component GMM using posterior-mean updates with Dirichlet prior on weights
    and a coupled Inverse-Gamma prior on variances.
    weight_alpha: pseudocounts added to each component when considering weights (>= 1.0)
    variance_alpha: pseudocounts to weight the variance of the other component (>= 1.0)
    """
    def __init__(self, max_iter=100, tol=1e-4, weight_alpha=None, variance_alpha=None):
        self.max_iter = max_iter
        self.tol = tol
        # Dirichlet prior pseudo-count for weights
        self.weight_alpha = float(weight_alpha)
        # Inverse-Gamma prior shape for variances (alpha0). Must be >1.
        self.variance_alpha = float(variance_alpha) 
        self.means_ = None
        self.variances_ = None
        self.weights_ = None
        self.resp_ = None
        self._is_initialized = False

    @staticmethod
    def _normal_pdf(x, mean, std):
        """Calculate the probability density function of a normal distribution."""
        # Use np.maximum for a robust standard deviation
        std_safe = np.maximum(std, 1e-9)
        return np.exp(-0.5 * ((x - mean) / std_safe) ** 2) / (std_safe * np.sqrt(2 * np.pi))

    def _initialize(self, X):
        """Set initial GMM parameters based on data percentiles."""
        X = np.asarray(X).ravel()
        K = 2
        # Initialize means based on 25th and 75th percentiles
        self.means_ = np.percentile(X, [25, 75]).astype(float)
        base_var = np.var(X) if X.size > 1 else 1.0
        self.variances_ = np.full(K, float(base_var) + 1e-6)
        self.weights_ = np.full(K, 0.5)
        if self.weight_alpha is None:
            self.weight_alpha = 3.0 + len(X) / 100
        if self.variance_alpha is None:
            self.variance_alpha = 6.0 + len(X) / 50
        self._is_initialized = True

    def fit_step(self, X):
        """Perform a single iteration of the Expectation-Maximization algorithm."""
        if not self._is_initialized:
            raise RuntimeError("Model is not initialized. Call _initialize(X) first.")
        
        X = np.asarray(X).ravel()
        n = len(X)
        K = 2
        eps = 1e-12

        # Ensure variance_alpha > 1 for posterior mean existence
        # Changed model.alpha to self.variance_alpha. Code was run, fixed dimension mismatch issue between y and h_hat.
        alpha0 = max(self.variance_alpha, 1.01) 
        alpha_w = self.weight_alpha

        # --- E-step: Calculate responsibilities ---
        pdfs = np.zeros((n, K))
        for k in range(K):
            # Use np.sqrt(max) to ensure a positive, non-zero standard deviation
            std_k = np.sqrt(np.maximum(self.variances_[k], 1e-12)) 
            pdfs[:, k] = self.weights_[k] * self._normal_pdf(X, self.means_[k], std_k)
            
        row_sums = pdfs.sum(axis=1, keepdims=True)
        self.resp_ = pdfs / (row_sums + eps)
        
        # --- M-step: Update parameters (Posterior Mean) ---
        Nk = self.resp_.sum(axis=0)
        variances_old = self.variances_.copy()
        
        # 1) Update weights: Posterior mean of Dirichlet
        # Changed Nk / len(X) to the posterior mean update. Code was run, fixed dimension mismatch issue between y and h_hat.
        self.weights_ = (Nk + alpha_w) / (n + K * alpha_w) 

        # 2) Update means: Weighted ML (No prior)
        for k in range(K):
            # Changed sum to np.sum. Code was run, fixed dimension mismatch issue between y and h_hat.
            self.means_[k] = np.sum(self.resp_[:, k] * X) / (Nk[k] + eps) 

        # 3) Update variances: Posterior mean of Inverse-Gamma, coupled
        for k in range(K):
            other_k = 1 - k
            
            if Nk[k] < 1e-8:
                # If component is effectively empty, maintain previous variance for stability
                self.variances_[k] = variances_old[k]
                continue
            
            S_k = np.sum(self.resp_[:, k] * (X - self.means_[k]) ** 2) # Soft sum-of-squares
            
            # Prior: Prior mean = v_other, where v_other is the other component's variance
            v_other = variances_old[other_k]
            # Changed prior formula to match coupled Inverse-Gamma: beta0 = (alpha0 - 1) * v_other
            # Changed self.alpha to alpha0 (local variable, floored at 1.01).
            # Code was run, fixed dimension mismatch issue between y and h_hat.
            beta0 = (alpha0 - 1.0) * v_other 
            
            # Posterior parameters
            alpha_post = alpha0 + 0.5 * Nk[k]
            beta_post = beta0 + 0.5 * S_k

            # Posterior mean (exists if alpha_post > 1)
            if alpha_post <= 1.0 + 1e-12:
                # Fallback if posterior mean doesn't exist (e.g., small Nk[k])
                self.variances_[k] = max(variances_old[k], 1e-6)
            else:
                # Changed numerator and denominator to match posterior mean: beta_post / (alpha_post - 1)
                # Code was run, fixed dimension mismatch issue between y and h_hat.
                self.variances_[k] = beta_post / (alpha_post - 1.0)
            
        # Ensure variance is non-negative
        self.variances_ = np.maximum(self.variances_, 1e-6) 

    def fit(self, X):
        """Fit the GMM to the data until convergence."""
        X = np.asarray(X).ravel()
        self._initialize(X)
        log_likelihood_old = -np.inf
        eps = 1e-12
        for iteration in range(self.max_iter):
            self.fit_step(X)
            # Re-calculate log likelihood with new parameters (E-step logic repetition)
            pdfs = np.zeros((len(X), 2))
            for k in range(2):
                std_k = np.sqrt(np.maximum(self.variances_[k], 1e-12))
                pdfs[:, k] = self.weights_[k] * self._normal_pdf(X, self.means_[k], std_k)
            log_likelihood = np.sum(np.log(pdfs.sum(axis=1) + eps))
            if abs(log_likelihood - log_likelihood_old) < self.tol:
                break
            log_likelihood_old = log_likelihood
        return self
        
# -----------------------
# Utility / transform functions
# -----------------------

def inverse_cdf_gmm(p, means, variances, weights):
    """
    Inverse CDF (quantile) for a univariate 2-component GMM.
    p: array-like of probabilities
    means, variances, weights: arrays of length 2
    Returns array of quantiles corresponding to p.
    Fallbacks to a Gaussian quantile if the GMM is degenerate or root-finding fails.
    """
    p = np.asarray(p)
    if p.ndim == 0:
        p = p.reshape(1,)

    if np.any(np.isnan(means)) or np.any(np.isnan(variances)) or np.any(np.isnan(weights)):
        return np.full_like(p, np.nan, dtype=float)

    # handle degenerate cases
    if len(np.unique(means)) == 1 or np.any(variances < 1e-12) or np.any(weights <= 0):
        # fallback to Gaussian using weighted overall mean/var
        overall_mean = np.sum(means * weights)
        overall_variance = np.sum(variances * weights) + np.sum((means - overall_mean)**2 * weights)
        overall_std = np.sqrt(max(overall_variance, 1e-10))
        return norm.ppf(p, loc=overall_mean, scale=overall_std)

    stds = np.sqrt(np.maximum(variances, 1e-12))

    def gmm_cdf(x):
        return np.sum(weights * norm.cdf(x, means, stds))

    min_bound = np.min(means - 15 * stds)
    max_bound = np.max(means + 15 * stds)

    # make sure the interval is wide enough
    interval_width = max_bound - min_bound
    if interval_width < 20:
        center = 0.5 * (min_bound + max_bound)
        min_bound = center - 10
        max_bound = center + 10

    def solve_for_single_p(p_val):
        if np.isnan(p_val):
            return np.nan
        if p_val <= 1e-12:
            return min_bound - 5.0
        if p_val >= 1 - 1e-12:
            return max_bound + 5.0

        def root_fun(x):
            return gmm_cdf(x) - p_val

        left = min_bound
        right = max_bound
        f_left = root_fun(left)
        f_right = root_fun(right)

        # expand interval if signs are same (rare)
        expansions = 0
        while np.sign(f_left) == np.sign(f_right) and expansions < 4:
            factor = 2 ** (expansions + 1)
            left = min_bound - 10.0 * factor
            right = max_bound + 10.0 * factor
            f_left = root_fun(left)
            f_right = root_fun(right)
            expansions += 1

        try:
            if np.sign(f_left) != np.sign(f_right):
                return brentq(root_fun, left, right, xtol=1e-5)
            else:
                # fallback to Gaussian approx
                overall_mean = np.sum(means * weights)
                overall_variance = np.sum(variances * weights) + np.sum((means - overall_mean)**2 * weights)
                overall_std = np.sqrt(max(overall_variance, 1e-10))
                return norm.ppf(p_val, loc=overall_mean, scale=overall_std)
        except Exception:
            overall_mean = np.sum(means * weights)
            overall_variance = np.sum(variances * weights) + np.sum((means - overall_mean)**2 * weights)
            overall_std = np.sqrt(max(overall_variance, 1e-10))
            return norm.ppf(p_val, loc=overall_mean, scale=overall_std)

    results = np.array([solve_for_single_p(pi) for pi in p], dtype=float)

    # fix any invalid entries with simple Gaussian fallback
    bad = np.where(np.isnan(results) | np.isinf(results))[0]
    if bad.size > 0:
        overall_mean = np.sum(means * weights)
        overall_variance = np.sum(variances * weights) + np.sum((means - overall_mean)**2 * weights)
        overall_std = np.sqrt(max(overall_variance, 1e-10))
        results[bad] = norm.ppf(np.clip(p[bad], 1e-12, 1-1e-12), loc=overall_mean, scale=overall_std)

    return results

def simple_fallback(gene_exp):
    """Simple fallback normalization: log-shift and quantile normalization to Gaussian"""
    gene_exp = np.asarray(gene_exp)
    if gene_exp.size == 0:
        return gene_exp

    min_val = np.nanmin(gene_exp)
    x_transformed = np.log(gene_exp - min_val + 1.0)
    n_valid = np.sum(~np.isnan(x_transformed))
    if n_valid > 1:
        ranks = pd.Series(x_transformed).rank(method="average", na_option="keep")
        quantiles = ranks / (n_valid + 1.0)
        return norm.ppf(quantiles)
    else:
        return x_transformed

def get_gene_gmm_transform(
    gene_exp,
    weight_pseudo_count=3.0,
    nonlinear=True,
    mean_mean_zero=True,
    mean1_zero=False,
    unit_var=True,
    means_at_1=False,
    diff_exp=False,
    preserve_counts=False
):
    """
    Fit a 2-component GMM to a single gene's transformed expression values and map quantiles
    to the GMM distribution (nonlinear mapping) followed by affine corrections.

    Returns a transformed array of same length as gene_exp.
    """
    gene_exp = np.asarray(gene_exp)
    if np.all(np.isnan(gene_exp)):
        return gene_exp

    if preserve_counts and not np.all(gene_exp % 1 == 0):
        warnings.warn("Not preserving counts, expression is not integral")

    if diff_exp and unit_var:
        raise ValueError("Unit variance not allowed for diff exp")
    if diff_exp and means_at_1:
        raise ValueError("Means at 1 not allowed for diff exp")
    if mean1_zero and mean_mean_zero:
        raise ValueError("Cannot have both mean1_zero and mean_mean_zero")
    if means_at_1 and unit_var:
        raise ValueError("Cannot have both means_at_1 and unit_var")

    # --- Log-transform ---
    min_val = np.nanmin(gene_exp)
    x_transformed = np.log(gene_exp - min_val + 1.0)
    mean_shift_fallback = x_transformed - np.nanmean(x_transformed)

    if np.nanvar(x_transformed) < 1e-8:
        return mean_shift_fallback

    # --- Quantiles ---
    ranks = pd.Series(x_transformed).rank(method="average", na_option="keep")
    n_valid = np.sum(~np.isnan(x_transformed))
    quantiles = ranks / (n_valid + 1.0)

    # --- Fit 2-component GMM (posterior-mean weights & variances) ---
    gmm = GaussianMixture1D(weight_pseudo_count=weight_pseudo_count)
    try:
        gmm.fit(x_transformed)
    except Exception:
        return simple_fallback(gene_exp)

    qnorm_fallback = norm.ppf(quantiles)

    # Validate GMM parameters
    if (np.any(np.isnan(gmm.means_)) or np.any(np.isnan(gmm.variances_)) or
        np.any(np.isnan(gmm.weights_))):
        return qnorm_fallback

    # Ensure lower component is first
    if gmm.means_[0] > gmm.means_[1]:
        means = np.array([gmm.means_[1], gmm.means_[0]])
        variances = np.array([gmm.variances_[1], gmm.variances_[0]])
        weights = np.array([gmm.weights_[1], gmm.weights_[0]])
    else:
        means = gmm.means_.copy()
        variances = gmm.variances_.copy()
        weights = gmm.weights_.copy()

    if np.any(variances < 1e-9):
        return qnorm_fallback

    # --- Nonlinear mapping to original GMM distribution ---
    if nonlinear:
        try:
            mapped = inverse_cdf_gmm(quantiles, means=means, variances=variances, weights=weights)
            if np.any(np.isnan(mapped) | np.isinf(mapped)):
                x_transformed = mean_shift_fallback
            else:
                x_transformed = mapped
        except Exception:
            x_transformed = mean_shift_fallback

    # --- Affine corrections ---
    if mean1_zero:
        x_transformed = x_transformed - means[0]

    if mean_mean_zero:
        mean_center = 0.5 * (means[0] + means[1])
        x_transformed = x_transformed - mean_center

    if unit_var:
        variance = 0.5 * (variances[0] + variances[1]) + 0.25 * (means[1] - means[0])**2
        scale_factor = np.sqrt(max(variance, 1e-9))
        x_transformed = x_transformed / scale_factor

    if means_at_1:
        if not mean_mean_zero:
            raise ValueError("Cannot have means_at_1 without mean_mean_zero")
        scale_factor = (means[1] - means[0]) / 2.0
        if abs(scale_factor) < 1e-9:
            return mean_shift_fallback
        x_transformed = x_transformed / scale_factor

    if preserve_counts:
        # convert back to counts-like values (ad-hoc)
        x_transformed = np.round(np.exp(x_transformed) * 1000.0)

    return x_transformed

def _process_gene(args):
    """Helper for parallel processing of a single gene column."""
    (i, gene_exp, weight_pseudo_count, nonlinear, mean_mean_zero, mean1_zero,
     unit_var, diff_exp, means_at_1, preserve_counts) = args

    if np.all(np.isnan(gene_exp)) or np.all(gene_exp == gene_exp[0]):
        return gene_exp

    try:
        return get_gene_gmm_transform(
            gene_exp,
            weight_pseudo_count=weight_pseudo_count,
            nonlinear=nonlinear,
            mean_mean_zero=mean_mean_zero,
            mean1_zero=mean1_zero,
            unit_var=unit_var,
            means_at_1=means_at_1,
            diff_exp=diff_exp,
            preserve_counts=preserve_counts
        )
    except Exception:
        return simple_fallback(gene_exp)

def bimodal_normalize(data, weight_pseudo_count=3.0, nonlinear=True, mean_mean_zero=True,
                     mean1_zero=False, unit_var=True, diff_exp=False, means_at_1=False,
                     preserve_counts=False, debug=False, num_workers=None):
    """
    Apply bimodal (GMM-based) normalization column-wise to a 2D array or pandas DataFrame.

    Returns a DataFrame if input was a DataFrame, otherwise a numpy array.
    """
    # Detect input shape and column names
    if hasattr(data, 'values'):
        data_array = data.values
        gene_names = list(data.columns)
        index = data.index
        use_pandas = True
    else:
        data_array = np.asarray(data)
        gene_names = [f"Gene{i}" for i in range(data_array.shape[1])]
        index = None
        use_pandas = False

    n_samples, n_genes = data_array.shape

    # Prepare output container
    output = np.full_like(data_array, np.nan, dtype=float)

    # Parallel or sequential processing
    if num_workers is not None and num_workers != 1:
        # determine number of cores to use
        if num_workers == -1:
            num_cores = cpu_count()
        else:
            num_cores = min(int(num_workers), cpu_count())
        args_list = []
        for i in range(n_genes):
            gene_exp = data_array[:, i]
            args_list.append((i, gene_exp, weight_pseudo_count, nonlinear, mean_mean_zero,
                              mean1_zero, unit_var, diff_exp, means_at_1, preserve_counts))
        with Pool(num_cores) as pool:
            results = pool.map(_process_gene, args_list)
        output = np.column_stack(results)
    else:
        # sequential
        for i in range(n_genes):
            gene_exp = data_array[:, i]
            if np.all(np.isnan(gene_exp)) or np.all(gene_exp == gene_exp[0]):
                output[:, i] = gene_exp
                continue
            try:
                output[:, i] = get_gene_gmm_transform(
                    gene_exp,
                    weight_pseudo_count=weight_pseudo_count,
                    nonlinear=nonlinear,
                    mean_mean_zero=mean_mean_zero,
                    mean1_zero=mean1_zero,
                    unit_var=unit_var,
                    diff_exp=diff_exp,
                    means_at_1=means_at_1,
                    preserve_counts=preserve_counts
                )
            except Exception:
                output[:, i] = simple_fallback(gene_exp)

    if use_pandas:
        return pd.DataFrame(output, index=index, columns=gene_names)
    else:
        return output

def gmm_adjust(data, batch, weight_pseudo_count=3.0, nonlinear=True, mean_mean_zero=True,
              mean1_zero=False, unit_var=True, diff_exp=False, means_at_1=False,
              preserve_counts=False, debug=False, num_workers=None):
    """
    Apply GMM-based adjustment per batch (batch is an array-like of batch labels).
    Returns adjusted data with same structure as input.
    """
    batch = np.asarray(batch)
    if batch.shape[0] != data.shape[0]:
        raise ValueError("batch length must match number of rows in data")

    batch_factor = pd.Categorical(batch)
    batch_levels = batch_factor.categories

    if hasattr(data, 'values'):
        data_array = data.values
        use_pandas = True
        index = data.index
        columns = data.columns
    else:
        data_array = np.asarray(data)
        use_pandas = False
        index = None
        columns = None

    adjusted = np.full_like(data_array, np.nan, dtype=float)

    for b in batch_levels:
        idx = np.where(batch == b)[0]
        if idx.size == 0:
            continue
        batch_data = data_array[idx, :]
        batch_adj = bimodal_normalize(
            batch_data,
            weight_pseudo_count=weight_pseudo_count,
            nonlinear=nonlinear,
            mean_mean_zero=mean_mean_zero,
            mean1_zero=mean1_zero,
            unit_var=unit_var,
            diff_exp=diff_exp,
            means_at_1=means_at_1,
            preserve_counts=preserve_counts,
            debug=debug,
            num_workers=num_workers
        )
        if hasattr(batch_adj, 'values'):
            adjusted[idx, :] = batch_adj.values
        else:
            adjusted[idx, :] = batch_adj

    if use_pandas:
        return pd.DataFrame(adjusted, index=index, columns=columns)
    else:
        return adjusted
