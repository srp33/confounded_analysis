import numpy as np
import pandas as pd
from scipy.stats import norm
from scipy.optimize import root_scalar
from tqdm import tqdm


class GaussianMixture1D:
    def __init__(self, n_components=2, max_iter=100, tol=1e-4, alpha0=1.0):
        self.n_components = n_components
        self.max_iter = max_iter
        self.tol = tol
        self.alpha0 = alpha0
        self.means_ = None
        self.variances_ = None
        self.weights_ = None
        self.resp_ = None
        self.log_likelihood_ = -np.inf

    def fit(self, X):
        X = X.ravel()
        n, K = X.shape[0], self.n_components
        eps = 1e-12

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

        self.means_, self.variances_, self.weights_ = means, variances, weights
        self.resp_, self.log_likelihood_ = responsibilities, log_likelihood
        return self

    def predict_proba(self, X):
        X = X.ravel()
        pdfs = np.array([
            self.weights_[k] * self._normal_pdf(X, self.means_[k], np.sqrt(self.variances_[k]))
            for k in range(self.n_components)
        ]).T
        return pdfs / (pdfs.sum(axis=1, keepdims=True) + 1e-12)

    @staticmethod
    def _normal_pdf(x, mean, std):
        return np.exp(-0.5 * ((x - mean) / std) ** 2) / (std * np.sqrt(2 * np.pi))


def inverse_cdf_gmm(p, means, variances, weights):
    p = np.asarray(p)
    eps = 1e-12

    if np.any(np.isnan(means)) or np.any(np.isnan(variances)) or np.any(np.isnan(weights)):
        return np.full_like(p, np.nan, dtype=float)

    if len(np.unique(means)) == 1 or np.any(variances < 1e-12):
        return norm.ppf(p)

    def gmm_cdf(x):
        return np.sum(weights * norm.cdf(x, loc=means, scale=np.sqrt(np.maximum(variances, 1e-10))))

    stds = np.sqrt(np.maximum(variances, 1e-10))
    min_bound, max_bound = np.min(means - 15 * stds), np.max(means + 15 * stds)

    if max_bound - min_bound < 20:
        center = 0.5 * (min_bound + max_bound)
        min_bound, max_bound = center - 10, center + 10

    result = np.empty_like(p)
    for i, pi in enumerate(p):
        if np.isnan(pi):
            result[i] = np.nan
            continue
        if pi <= 1e-10:
            result[i] = min_bound - 5
            continue
        if pi >= 1 - 1e-10:
            result[i] = max_bound + 5
            continue

        def f(x):
            return gmm_cdf(x) - pi

        try:
            sol = root_scalar(f, bracket=[min_bound, max_bound], method="brentq", xtol=1e-4)
            result[i] = sol.root if sol.converged else norm.ppf(pi)
        except Exception:
            result[i] = norm.ppf(pi)

    result[np.isinf(result) | np.isnan(result)] = norm.ppf(p[np.isinf(result) | np.isnan(result)])
    return result


def get_gene_gmm_transform(
    gene_exp,
    alpha0=10,
    nonlinear=True,
    mean_mean_zero=True,
    unit_var=True,
    means_at_1=False
):
    gene_exp = np.asarray(gene_exp, dtype=float)
    if np.all(np.isnan(gene_exp)):
        return gene_exp

    min_val = np.nanmin(gene_exp)
    x_transformed = np.log(gene_exp - min_val + 1)
    if np.nanvar(x_transformed) < 1e-8:
        return x_transformed - np.nanmean(x_transformed)

    valid_mask = ~np.isnan(x_transformed)
    ranks = np.argsort(np.argsort(x_transformed[valid_mask])) + 1
    quantiles = np.empty_like(x_transformed)
    quantiles[:] = np.nan
    quantiles[valid_mask] = ranks / (np.sum(valid_mask) + 1)

    gmm = GaussianMixture1D(n_components=2, alpha0=alpha0)
    gmm.fit(x_transformed[valid_mask])

    means, variances, weights = gmm.means_, gmm.variances_, gmm.weights_

    if means[0] > means[1]:
        means, variances, weights = means[::-1], variances[::-1], weights[::-1]

    if np.any(np.isnan(means)) or np.any(variances < 1e-9):
        mapped = norm.ppf(quantiles)
    elif nonlinear:
        mapped = inverse_cdf_gmm(quantiles, means, variances, weights)
    else:
        mapped = x_transformed

    if mean_mean_zero:
        mapped -= 0.5 * np.sum(means)

    if unit_var:
        variance = 0.5 * np.sum(variances) + 0.25 * (means[1] - means[0]) ** 2
        scale_factor = np.sqrt(max(variance, 1e-9))
        mapped /= scale_factor

    if means_at_1:
        if unit_var or not mean_mean_zero:
            raise ValueError("Cannot combine means_at_1 with unit_var or mean_mean_zero=False")
        scale_factor = (means[1] - means[0]) / 2
        mapped /= max(scale_factor, 1e-6)

    return mapped


def bimodal_normalize(data, alpha0=10, mean_only=False, debug=False):
    genes = data.columns
    transformed = pd.DataFrame(index=data.index, columns=genes, dtype=float)

    for i, gene in enumerate(genes, start=1):
        if debug and i % 1000 == 0:
            print(f"Processed {i}/{len(genes)} genes")

        gene_exp = data[gene].values
        try:
            transformed[gene] = get_gene_gmm_transform(
                gene_exp, alpha0=alpha0, nonlinear=not mean_only
            )
        except Exception as e:
            if debug:
                print(f"Gene {gene}: fallback due to {e}")
            valid = ~np.isnan(gene_exp)
            ranks = np.argsort(np.argsort(gene_exp[valid])) + 1
            q = ranks / (np.sum(valid) + 1)
            fallback = norm.ppf(q)
            result = np.full_like(gene_exp, np.nan)
            result[valid] = fallback / np.nanstd(fallback)
            transformed[gene] = result

    return transformed


def gmm_adjust(data, batch, alpha0=10, mean_only=False, debug=False):
    if debug:
        print(f"Starting GMM adjustment with {len(np.unique(batch))} batches")

    adjusted = pd.DataFrame(index=data.index, columns=data.columns, dtype=float)
    for b in np.unique(batch):
        if debug:
            print(f"Processing batch: {b}")
        idx = np.where(batch == b)[0]
        batch_data = data.iloc[idx]
        adjusted.iloc[idx] = bimodal_normalize(batch_data, alpha0, mean_only, debug)

    return adjusted