import numpy as np
import pandas as pd
from scipy.stats import norm
from tqdm import tqdm

class GaussianMixture1D:
    def __init__(self, n_components=2, max_iter=100, tol=1e-4, alpha0=10.0):
        self.n_components = n_components
        self.max_iter = max_iter
        self.tol = tol
        self.alpha0 = alpha0
        self.means_ = None
        self.variances_ = None
        self.weights_ = None
        self.resp_ = None

    @staticmethod
    def _normal_pdf(x, mean, std):
        return np.exp(-0.5 * ((x - mean) / std) ** 2) / (std * np.sqrt(2 * np.pi))

    def fit(self, X):
        X = np.asarray(X).ravel()
        n = len(X)
        K = self.n_components
        eps = 1e-12

        # Initialize means using percentiles
        percentiles = np.linspace(0, 100, K + 2)[1:-1]
        means = np.percentile(X, percentiles)
        variances = np.full(K, np.var(X))
        weights = np.full(K, 1.0 / K)
        log_likelihood_old = -np.inf

        for _ in range(self.max_iter):
            pdfs = np.array([
                weights[k] * self._normal_pdf(X, means[k], np.sqrt(variances[k]))
                for k in range(K)
            ]).T
            responsibilities = pdfs / (pdfs.sum(axis=1, keepdims=True) + eps)
            Nk = responsibilities.sum(axis=0)

            # MAP update for weights
            weights = (Nk + self.alpha0 - 1) / (n + K * (self.alpha0 - 1))
            means = (responsibilities.T @ X) / (Nk + eps)
            variances = np.array([
                np.sum(responsibilities[:, k] * (X - means[k]) ** 2) / (Nk[k] + eps)
                for k in range(K)
            ])
            variances = np.maximum(variances, 1e-6)

            log_likelihood = np.sum(np.log(pdfs.sum(axis=1) + eps))
            if abs(log_likelihood - log_likelihood_old) < self.tol:
                break
            log_likelihood_old = log_likelihood

        self.means_ = means
        self.variances_ = variances
        self.weights_ = weights
        self.resp_ = responsibilities
        return self

    def predict_proba(self, X):
        X = np.asarray(X).ravel()
        pdfs = np.array([
            self.weights_[k] * self._normal_pdf(X, self.means_[k], np.sqrt(self.variances_[k]))
            for k in range(self.n_components)
        ]).T
        return pdfs / (pdfs.sum(axis=1, keepdims=True) + 1e-12)

def inverse_cdf_gmm(p, means, variances, weights):
    def gmm_cdf(x):
        return np.sum(weights * norm.cdf(x, means, np.sqrt(variances)))

    results = []
    for p_val in p:
        if p_val <= 1e-10:
            results.append(min(means) - 5)
            continue
        if p_val >= 1 - 1e-10:
            results.append(max(means) + 5)
            continue
        try:
            from scipy.optimize import brentq
            lower = min(means) - 10 * np.sqrt(max(variances))
            upper = max(means) + 10 * np.sqrt(max(variances))
            root = brentq(lambda x: gmm_cdf(x) - p_val, lower, upper)
            results.append(root)
        except Exception:
            results.append(norm.ppf(p_val))
    return np.array(results)

def get_gene_gmm_transform(
    gene_exp,
    alpha0=10,
    nonlinear=True,
    mean_mean_zero=True,
    mean1_zero=False,
    unit_var=True,
    means_at_1=False,
    diff_exp=False,
    preserve_counts=False
):
    gene_exp = np.asarray(gene_exp)
    if np.all(np.isnan(gene_exp)):
        return gene_exp

    if preserve_counts and not np.all(np.floor(gene_exp) == gene_exp):
        print("Warning: Not preserving counts, expression not integral")

    if diff_exp and unit_var:
        raise ValueError("Unit variance not allowed for diff exp")
    if diff_exp and means_at_1:
        raise ValueError("Means at 1 not allowed for diff exp")

    # Log transform
    min_val = np.nanmin(gene_exp)
    x_transformed = np.log(gene_exp - min_val + 1)
    mean_shift_fallback = x_transformed - np.nanmean(x_transformed)

    if np.nanvar(x_transformed) < 1e-8:
        return mean_shift_fallback

    # Quantiles
    ranks = pd.Series(x_transformed).rank(method="average", na_option="keep")
    n_valid = np.sum(~np.isnan(x_transformed))
    quantiles = ranks / (n_valid + 1)

    # Fit 2-component GMM
    gmm = GaussianMixture1D(alpha0=alpha0)
    gmm.fit(x_transformed)

    means = np.sort(gmm.means_)
    variances = gmm.variances_[np.argsort(gmm.means_)]
    weights = gmm.weights_[np.argsort(gmm.means_)]

    if np.any(np.isnan(means)) or np.any(variances < 1e-9):
        return norm.ppf(quantiles)

    if nonlinear:
        try:
            mapped = inverse_cdf_gmm(quantiles, means, variances, weights)
            if np.any(np.isnan(mapped)) or np.any(np.isinf(mapped)):
                x_transformed = mean_shift_fallback
            else:
                x_transformed = mapped
        except Exception:
            x_transformed = mean_shift_fallback

    # Affine corrections
    if diff_exp:
        if mean_mean_zero or unit_var or means_at_1:
            print("Warning: diff_exp disables centering and scaling")
        x_transformed = x_transformed - means[0]
        mean_mean_zero = unit_var = means_at_1 = False

    if mean1_zero:
        # Adjusts the first mean to be 0
        if mean_mean_zero:
            raise ValueError("Cannot have both mean1_zero and mean_mean_zero")
        if means_at_1:
            raise ValueError("Cannot have both mean1_zero and means_at_1")
        x_transformed = x_transformed - means[0]

    if mean_mean_zero:
        mean_center = 0.5 * (means[0] + means[1])
        x_transformed = x_transformed - mean_center

    if unit_var:
        variance = 0.5 * (variances[0] + variances[1]) + 0.25 * (means[1] - means[0]) ** 2
        scale_factor = np.sqrt(max(variance, 1e-9))
        x_transformed = x_transformed / scale_factor

    if means_at_1:
        if unit_var:
            raise ValueError("Cannot have both means_at_1 and unit_var")
        if not mean_mean_zero:
            raise ValueError("Cannot have means_at_1 without mean_mean_zero")
        scale_factor = (means[1] - means[0]) / 2
        if abs(scale_factor) < 1e-6:
            return mean_shift_fallback
        x_transformed = x_transformed / scale_factor

    if preserve_counts:
        exp_transformed = np.exp(x_transformed)
        x_transformed = np.round(exp_transformed * 1000)

    return x_transformed

def bimodal_normalize(data, alpha0=10, mean_only=False, diff_exp=False, preserve_counts=False, debug=False):
    genes = data.columns
    n_genes = len(genes)
    n_samples = data.shape[0]
    normalized = np.full_like(data, np.nan, dtype=float)

    for i, gene in enumerate(genes):
        if debug and i % 500 == 0:
            print(f"Processing {i}/{n_genes} genes")

        gene_exp = data.iloc[:, i].values

        if np.all(np.isnan(gene_exp)) or np.nanvar(gene_exp) < 1e-8:
            normalized[:, i] = gene_exp
            continue

        try:
            normed = get_gene_gmm_transform(
                gene_exp,
                alpha0=alpha0,
                nonlinear=not mean_only,
                diff_exp=diff_exp,
                preserve_counts=preserve_counts
            )
            normalized[:, i] = normed
        except Exception as e:
            if debug:
                print(f"Error on {gene}: {e}")
            normalized[:, i] = gene_exp

    return pd.DataFrame(normalized, index=data.index, columns=data.columns)

def gmm_adjust(data, batch, alpha0=10, mean_only=False, diff_exp=False, preserve_counts=False, debug=False):
    if debug:
        print(f"Starting GMM adjustment with {len(np.unique(batch))} batches")

    batch = pd.Series(batch, index=data.index)
    adjusted = np.full_like(data, np.nan, dtype=float)

    for b in batch.unique():
        if debug:
            print(f"Processing batch {b}")
        batch_mask = batch == b
        batch_data = data.loc[batch_mask]
        adjusted_batch = bimodal_normalize(
            batch_data,
            alpha0=alpha0,
            mean_only=mean_only,
            diff_exp=diff_exp,
            preserve_counts=preserve_counts,
            debug=debug
        )
        adjusted.loc[batch_mask] = adjusted_batch.values

    return pd.DataFrame(adjusted, index=data.index, columns=data.columns)
