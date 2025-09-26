import numpy as np

class GaussianMixture1D:
    def __init__(self, n_components=2, max_iter=100, tol=1e-4, alpha0=1.0):
        """
            n_components (int): The number of mixture components.
            max_iter (int): The maximum number of EM iterations.
            tol (float): The convergence tolerance for the log-likelihood.
            alpha0 (float): The concentration parameter for the symmetric Dirichlet prior on the mixture weights.
                alpha0 = 1.0 corresponds to no prior (MLE).
                alpha0 > 1.0 encourages weights to be similar (MAP).
        """
        self.n_components = n_components
        self.max_iter = max_iter
        self.tol = tol
        self.alpha0 = alpha0
        self.means_ = None
        self.variances_ = None
        self.weights_ = None
        self.log_likelihood_ = -np.inf

    def fit(self, X):
        X = X.ravel()
        n = X.shape[0]
        K = self.n_components
        eps = 1e-12

        # Initialize means using percentiles for better starting points
        if K == 1:
            means = np.array([X.mean()])
        else:
            # Spacing K points between 0 and 100
            percentiles = np.linspace(0, 100, K + 2)[1:-1]
            means = np.percentile(X, percentiles)

        variances = np.full(K, X.var())
        weights = np.full(K, 1.0 / K)

        log_likelihood_old = -np.inf

        for _ in range(self.max_iter):
            # E-step: responsibilities
            pdfs = np.array([
                weights[k] * self._normal_pdf(X, means[k], np.sqrt(variances[k]))
                for k in range(K)
            ])
            pdfs = pdfs.T  # shape (n, k)
            responsibilities = pdfs / (pdfs.sum(axis=1, keepdims=True) + eps)

            # M-step: update parameters
            Nk = responsibilities.sum(axis=0)

            # Update weights with Dirichlet prior (MAP estimate)
            # The original MLE update was: weights = Nk / n
            weights = (Nk + self.alpha0 - 1) / (n + K * (self.alpha0 - 1))

            means = (responsibilities.T @ X) / (Nk + eps)

            variances = np.array([
                np.sum(responsibilities[:, k] * (X - means[k]) ** 2) / (Nk[k] + 1e-12)
                for k in range(K)
            ])
            variances = np.maximum(variances, 1e-6) # Variance floor


            # Check convergence
            log_likelihood = np.sum(np.log(pdfs.sum(axis=1) + eps))
            if abs(log_likelihood - log_likelihood_old) < self.tol:
                break
            log_likelihood_old = log_likelihood

        self.means_ = means
        self.variances_ = variances
        self.weights_ = weights
        self.resp_ = responsibilities
        self.log_likelihood_ = log_likelihood
        return self

    def predict_proba(self, X):
        X = X.ravel()
        pdfs = np.array([
            self.weights_[k] * self._normal_pdf(X, self.means_[k], np.sqrt(self.variances_[k]))
            for k in range(self.n_components)
        ])
        pdfs = pdfs.T
        return pdfs / (pdfs.sum(axis=1, keepdims=True) + 1e-12)

    def bic(self, X):
        X = X.ravel()
        n, k = X.shape[0], self.n_components
        log_likelihood = self.log_likelihood_
        # Parameters: means (k), variances (k), weights (k-1)
        n_params = 2 * k + (k - 1)
        return -2 * log_likelihood + n_params * np.log(n)

    @staticmethod
    def _normal_pdf(x, mean, std):
        return np.exp(-0.5 * ((x - mean) / std) ** 2) / (std * np.sqrt(2 * np.pi))


class Gaussian1Component(GaussianMixture1D):
    def __init__(self, n_components=1):
        super().__init__(n_components=n_components, max_iter=0, tol=0)
        self.means_ = None
        self.variances_ = None
        self.weights_ = None
        self.resp_ = None
        self.log_likelihood_ = None
    
    def fit(self, X):
        # Since we only fit one component, this has a closed-form solution
        X = X.ravel()
        self.means_ = [X.mean()]
        self.variances_ = [X.var()]
        self.weights_ = [1.0]
        self.resp_ = np.ones((X.shape[0], 1))
        self.log_likelihood_ = np.sum(np.log(self._normal_pdf(X, self.means_[0], np.sqrt(self.variances_[0]))))
        return self
    
    def predict_proba(self, X):
        X = X.ravel()
        return np.ones((X.shape[0], 1))

    def bic(self, X):
        X = X.ravel()
        n = X.shape[0]
        log_likelihood = self.log_likelihood_
        # Parameters: means (1), variances (1)
        n_params = 2
        return -2 * log_likelihood + n_params * np.log(n)


    

def get_gmm_responsibilites(log_exp):
    column_names = log_exp.columns
    log_exp_np = log_exp.to_numpy()
    responsibilities_dict = {}
    
    for gene_exp, gene_name in tqdm(zip(log_exp_np.T, column_names)):
        # gene_exp = gene_exp.reshape(-1, 1)

        gmm = GaussianMixture1D(n_components=2, alpha0=10)
        gmm.fit(gene_exp)
        gmm1 = Gaussian1Component()
        gmm1.fit(gene_exp)

        # Skip unimodal genes
        # if gmm.bic(gene_exp) >= gmm1.bic(gene_exp):
        #     continue

        # Calulate the mixture responsibilites for each sample
        # These are later used to weight the samples
        # Only get the responsibilites for the lower component, to save space
        responsibilites = gmm.predict_proba(gene_exp)[:, 0] 

        # Make sure responsibilites are for the lower component
        if gmm.means_[0] > gmm.means_[1]:
            responsibilites = 1 - responsibilites

        responsibilities_dict[gene_name] = responsibilites
    
    return pd.DataFrame(responsibilities_dict)

responsibilites_20194 = get_gmm_responsibilites(log_exp_20194)