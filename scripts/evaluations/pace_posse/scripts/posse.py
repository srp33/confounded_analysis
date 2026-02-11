import numpy as np
from scipy.special import softmax
from scipy.stats import gamma
from dataclasses import dataclass
from typing import List, Dict, Tuple, Optional, Union
import os
import sys

import gseapy as gp
from sklearn.decomposition import PCA

# ==========================================
# Global Constants
# ==========================================
EPS = 1e-8

# ==========================================
# Configuration & Data Structures
# ==========================================

@dataclass
class POSSEHyperparameters:
    """
    Hyperparameters for POSSE v5.0 (ComBat-Initialized).
    Uses ComBat estimates as gene-specific priors for biological refinement.
    """
    tau: float = 25.0          # Sharper contrast than v4.2
    top_k_percent: float = 0.20   # Broader context to stabilize local estimates
    eta: float = 0.3           # Faster adaptation
    max_iter: int = 5          # Iterations to refine the ComBat prior
    min_pathway_size: int = 5  # Minimum genes required to form a pathway context
    hk_percentile: float = 0.2 # Housekeeping gene percentile (bottom 20% CoV)
    # Two-phase parameters (optional)
    phase1_tau: float = None   # Conservative tau for prior learning (if None, use tau)
    phase1_top_k: float = None # Conservative top-K for prior learning (if None, use top_k_percent)
    

@dataclass
class BatchData:
    """
    Container for batch expression data.
    data: Matrix of shape (Genes, Samples)
    gene_indices: Array of gene symbols/IDs corresponding to rows
    """
    data: np.ndarray 
    gene_indices: np.ndarray

# ==========================================
# Math Utilities
# ==========================================

def centered_cosine_similarity(X: np.ndarray, Y: np.ndarray) -> np.ndarray:
    """
    Compute centered cosine similarity (Pearson correlation) between samples.
    Used for Shape-based Navigation.
    
    Args:
        X: Reference data (Genes, N)
        Y: Target data (Genes, M)
    Returns: 
        Similarity matrix (N, M)
    """
    # Center data gene-wise
    X_mean = np.mean(X, axis=0, keepdims=True)
    Y_mean = np.mean(Y, axis=0, keepdims=True)
    
    X_c = X - X_mean
    Y_c = Y - Y_mean
    
    # Compute norms
    X_norm = np.linalg.norm(X_c, axis=0, keepdims=True) + EPS
    Y_norm = np.linalg.norm(Y_c, axis=0, keepdims=True) + EPS
    
    # Cosine calculation
    return (X_c / X_norm).T @ (Y_c / Y_norm)

def safe_entropy(p):
    """Entropy of columns (Target samples)"""
    return -np.sum(p * np.log(p + EPS), axis=0)

def weighted_stats(data: np.ndarray, weights: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute weighted mean and variance along axis 1 (samples).
    
    Args:
        data: (Genes, Samples)
        weights: (Samples,) - Normalized weights summing to ~1
    Returns:
        mu: (Genes,)
        var: (Genes,)
    """
    w_sum = np.sum(weights) + EPS
    
    # Weighted Mean
    mu = np.sum(data * weights, axis=1) / w_sum
    
    # Weighted Variance: sum(w * (x - mu)^2) / sum(w)
    centered = data - mu[:, np.newaxis]
    var = np.sum(weights * (centered**2), axis=1) / w_sum
    
    return mu, var

# ==========================================
# Legacy Class for Compatibility
# ==========================================
class ComBatBaseline:
    """ComBat baseline using pure Python implementation for gene-specific prior initialization in POSSE v5.0"""
    
    def compute_baseline(self, X, Y):
        """
        Compute gene-specific ComBat priors using pure Python ComBat implementation.
        
        Args:
            X: Reference data (genes x samples)
            Y: Target data (genes x samples)
            
        Returns:
            alpha: Per-gene scale factors (genes,)
            beta: Per-gene shift factors (genes,)
        """
        return self._compute_combat_python(X, Y)
    
    def _compute_combat_python(self, X, Y):
        """
        Pure Python implementation of ComBat algorithm based on Johnson et al. 2007.
        """
        from scipy.stats import gamma
        from scipy.optimize import minimize_scalar
        
        # DIAGNOSTIC: Print data characteristics before ComBat
        print(f"[COMBAT DIAGNOSTIC] Input data characteristics:")
        print(f"  X (reference) shape: {X.shape}")
        print(f"  Y (target) shape: {Y.shape}")
        print(f"  X range: [{np.min(X):.3f}, {np.max(X):.3f}]")
        print(f"  Y range: [{np.min(Y):.3f}, {np.max(Y):.3f}]")
        print(f"  X mean: {np.mean(X):.3f}, std: {np.std(X):.3f}")
        print(f"  Y mean: {np.mean(Y):.3f}, std: {np.std(Y):.3f}")
        
        # Combine data: X (reference, batch=1) and Y (target, batch=2)
        combined_data = np.hstack([X, Y])  # Shape: (genes, X_samples + Y_samples)
        batch_vector = np.concatenate([
            np.ones(X.shape[1], dtype=int),      # Reference batch = 1
            np.full(Y.shape[1], 2, dtype=int)    # Target batch = 2
        ])
        
        n_genes, n_samples = combined_data.shape
        n_batches = 2
        batches = [np.where(batch_vector == 1)[0], np.where(batch_vector == 2)[0]]
        n_batch_samples = [len(batch) for batch in batches]
        
        print(f"  Combined data shape: {combined_data.shape}")
        print(f"  Batch sizes: {n_batch_samples}")
        
        # Filter zero-variance genes
        gene_vars = np.var(combined_data, axis=1)
        varying_genes = gene_vars > 1e-8
        
        if np.sum(varying_genes) < n_genes:
            print(f"  Filtering {n_genes - np.sum(varying_genes)} zero-variance genes for ComBat")
        
        # Work only with varying genes
        data_filt = combined_data[varying_genes, :]
        n_genes_filt = data_filt.shape[0]
        
        # Create design matrix (batch indicators without intercept)
        # For 2 batches, we need only 1 indicator (reference coding)
        design = np.zeros((n_samples, 1))
        design[batches[1], 0] = 1  # Target batch indicator (batch 2 = 1, batch 1 = 0)
        
        # Estimate coefficients B_hat using least squares
        # This gives us the difference between batches
        B_hat = np.linalg.solve(design.T @ design + 1e-8 * np.eye(design.shape[1]), 
                               design.T @ data_filt.T).T
        
        # Calculate grand mean (reference batch mean)
        grand_mean = np.mean(data_filt[:, batches[0]], axis=1)  # Reference batch mean
        
        # Calculate pooled variance using reference batch
        ref_residuals = data_filt[:, batches[0]] - grand_mean[:, np.newaxis]
        var_pooled = np.var(ref_residuals, axis=1, ddof=1)
        var_pooled = np.maximum(var_pooled, 1e-8)  # Avoid division by zero
        
        # Standardize data
        stand_mean = np.tile(grand_mean[:, np.newaxis], (1, n_samples))
        s_data = (data_filt - stand_mean) / np.sqrt(var_pooled)[:, np.newaxis]
        
        print(f"[COMBAT DIAGNOSTIC] Standardization completed")
        print(f"  Standardized data range: [{np.min(s_data):.3f}, {np.max(s_data):.3f}]")
        
        # Estimate batch effects (difference from reference)
        gamma_hat = np.zeros((n_genes_filt, n_batches))
        gamma_hat[:, 0] = 0  # Reference batch effect = 0
        gamma_hat[:, 1] = B_hat[:, 0]  # Target batch effect = difference from reference
        
        # Estimate variance parameters for each batch
        delta_hat = np.zeros((n_batches, n_genes_filt))
        for i, batch_idx in enumerate(batches):
            batch_residuals = s_data[:, batch_idx] - gamma_hat[:, i:i+1]
            delta_hat[i, :] = np.var(batch_residuals, axis=1, ddof=1)
            delta_hat[i, :] = np.maximum(delta_hat[i, :], 1e-8)
        
        # Empirical Bayes priors
        gamma_bar = np.mean(gamma_hat, axis=1)
        t2 = np.var(gamma_hat, axis=1, ddof=1)
        t2 = np.maximum(t2, 1e-8)
        
        # Estimate inverse gamma priors for delta
        def estimate_inv_gamma_params(x):
            """Estimate inverse gamma parameters using method of moments"""
            x = x[x > 1e-8]  # Remove near-zero values
            if len(x) < 2:
                return 1.0, 1.0
            
            mean_x = np.mean(x)
            var_x = np.var(x, ddof=1)
            
            # Method of moments for inverse gamma
            if var_x > 0 and mean_x > 0:
                # For inverse gamma: mean = b/(a-1), var = b^2/((a-1)^2(a-2))
                # Solving: a = 2 + mean^2/var, b = mean*(a-1)
                a = 2 + mean_x**2 / var_x
                b = mean_x * (a - 1)
                return max(a, 1.1), max(b, 1e-8)
            else:
                return 1.0, 1.0
        
        a_prior = np.zeros(n_batches)
        b_prior = np.zeros(n_batches)
        for i in range(n_batches):
            a_prior[i], b_prior[i] = estimate_inv_gamma_params(delta_hat[i, :])
        
        print(f"[COMBAT DIAGNOSTIC] Prior estimation completed")
        print(f"  Gamma priors - mean: {np.mean(gamma_bar):.3f}, var: {np.mean(t2):.3f}")
        print(f"  Delta priors - a: {a_prior}, b: {b_prior}")
        
        # Empirical Bayes estimation
        gamma_star = np.zeros_like(gamma_hat)
        delta_star = np.zeros_like(delta_hat)
        
        for i in range(n_batches):
            # Posterior mean for gamma (normal prior)
            precision = 1.0 / t2 + len(batches[i])
            gamma_star[:, i] = (gamma_bar / t2 + len(batches[i]) * gamma_hat[:, i]) / precision
            
            # Posterior parameters for delta (inverse gamma prior)
            for g in range(n_genes_filt):
                n_i = len(batches[i])
                ss = np.sum((s_data[g, batches[i]] - gamma_star[g, i])**2)
                
                a_post = a_prior[i] + n_i / 2.0
                b_post = b_prior[i] + ss / 2.0
                
                # Posterior mean of inverse gamma
                if a_post > 1:
                    delta_star[i, g] = b_post / (a_post - 1)
                else:
                    delta_star[i, g] = delta_hat[i, g]
        
        # Set reference batch parameters to zero/one
        gamma_star[:, 0] = 0  # Reference batch mean adjustment = 0
        delta_star[0, :] = 1  # Reference batch variance adjustment = 1
        
        print(f"[COMBAT DIAGNOSTIC] Empirical Bayes estimation completed")
        print(f"  Gamma* range: [{np.min(gamma_star):.3f}, {np.max(gamma_star):.3f}]")
        print(f"  Delta* range: [{np.min(delta_star):.3f}, {np.max(delta_star):.3f}]")
        
        # Extract transformation parameters directly from ComBat estimates
        # For target batch (batch 1, index 1), the correction is:
        # Y_corrected = (Y_standardized - gamma_star) / sqrt(delta_star)
        # Converting back: Y_corrected = alpha * Y + beta where:
        # alpha = 1 / sqrt(delta_star) * sqrt(var_pooled) / sqrt(var_pooled) = 1 / sqrt(delta_star)  
        # beta = -gamma_star / sqrt(delta_star) * sqrt(var_pooled) + grand_mean - grand_mean / sqrt(delta_star)
        
        # Initialize for all genes
        alpha = np.ones(X.shape[0])
        beta = np.zeros(X.shape[0])
        
        # Apply to varying genes only
        if np.sum(varying_genes) > 0:
            # Target batch is batch 1 (index 1)
            target_gamma = gamma_star[varying_genes, 1]  # Mean adjustment for target batch
            target_delta = delta_star[1, :]              # Variance adjustment for target batch
            target_var_pooled = var_pooled               # Pooled variance
            target_grand_mean = grand_mean               # Grand mean
            
            # Direct transformation parameters from ComBat correction
            # The correction transforms: (Y - grand_mean)/sqrt(var_pooled) -> ((Y - grand_mean)/sqrt(var_pooled) - gamma)/sqrt(delta)
            # Rearranging: Y_corrected = Y/sqrt(delta) - gamma*sqrt(var_pooled)/sqrt(delta) + grand_mean*(1 - 1/sqrt(delta))
            alpha[varying_genes] = 1.0 / np.sqrt(target_delta)
            beta[varying_genes] = (-target_gamma * np.sqrt(target_var_pooled) / np.sqrt(target_delta) + 
                                  target_grand_mean * (1.0 - 1.0 / np.sqrt(target_delta)))
        
        print(f"[COMBAT DIAGNOSTIC] Direct transformation parameters:")
        print(f"  Alpha range: [{np.min(alpha):.6f}, {np.max(alpha):.6f}]")
        print(f"  Beta range: [{np.min(beta):.3f}, {np.max(beta):.3f}]")
        print(f"  Alpha mean: {np.mean(alpha):.6f}, std: {np.std(alpha):.6f}")
        print(f"  Beta mean: {np.mean(beta):.3f}, std: {np.std(beta):.3f}")
        
        return alpha, beta

# ==========================================
# Main Algorithm: POSSE v4.0
# ==========================================

class POSSE:
    """
    POSSE v5.0: ComBat-Initialized Local Refinement.
    Uses actual ComBat estimates from sva library as gene-specific priors for biological refinement.
    
    Mechanism:
    1. ComBat Initialization: Calculate gene-specific priors using actual ComBat from sva library
    2. Local Refinement: Use pathway-based peer finding to refine or override priors
    3. Trust Gating: High trust preserves biology, low trust defaults to ComBat priors
    """
    
    def __init__(self, 
                 pathway_dict: Dict[str, List[str]] = None,
                 pathway_source: str = 'MSigDB_Hallmark_2020',
                 organism: str = 'Human',
                 hyperparams: POSSEHyperparameters = None,
                 debug: bool = False):
        
        self.hp = hyperparams or POSSEHyperparameters()
        self.debug = debug
        self.gene_stability_scores = None
        self.combat = ComBatBaseline()  # For prior initialization
        
        if pathway_dict:
            self.pathway_dict = pathway_dict
        else:
            self._load_pathways(pathway_source, organism)

    def _load_pathways(self, name: str, organism: str):
        # First try to load from cache - fix path calculation
        cache_file = os.path.join(os.path.dirname(__file__), 'pathways_cache.pkl')
        
        if self.debug: print(f"DEBUG: Looking for cached pathways at: {cache_file}")
        
        if os.path.exists(cache_file):
            try:
                import pickle
                if self.debug: print(f"DEBUG: Loading cached pathways from {cache_file}")
                with open(cache_file, 'rb') as f:
                    self.pathway_dict = pickle.load(f)
                if self.debug: print(f"DEBUG: Loaded {len(self.pathway_dict)} cached pathways")
                
                # FAIL FAST: Ensure pathways were actually loaded
                if not self.pathway_dict:
                    raise RuntimeError("CRITICAL: Cached pathway file is empty")
                return
            except Exception as e:
                if self.debug: print(f"DEBUG: Failed to load cached pathways: {e}")
        
        # Fallback to online download
        try:
            if self.debug: print(f"DEBUG: Downloading pathway set: {name}")
            self.pathway_dict = gp.get_library(name=name, organism=organism)
            
            # FAIL FAST: Ensure pathways were downloaded
            if not self.pathway_dict:
                raise RuntimeError("CRITICAL: No pathways downloaded from online source")
                
        except Exception as e:
            # FAIL FAST: Don't continue without pathways
            raise RuntimeError(f"CRITICAL: Cannot load pathways. {e}. "
                             f"Run 'python download_pathways.py' first to cache pathways locally, "
                             f"or ensure network connectivity for online download.")

    def calculate_gene_stability(self, X):
        """
        Identify genes that are stable (Housekeeping-like) in the Reference.
        Metric: Inverse Coefficient of Variation.
        Returns: Score 0.0 (Volatile) to 1.0 (Stable Anchor).
        """
        mu = np.mean(X, axis=1)
        sigma = np.std(X, axis=1)
        
        # CoV = Sigma / Mu
        # Filter low expression noise
        mask = mu > np.median(mu)
        
        cov = np.ones_like(mu) * 100.0 # Default high variance
        cov[mask] = sigma[mask] / (np.abs(mu[mask]) + EPS)
        
        # Rank: Lower CoV = Higher Rank
        # argsort twice gives rank (0 to N-1)
        ranks = np.argsort(np.argsort(cov))
        
        # Normalize to 0-1 (1.0 = Lowest CoV / Most Stable)
        stability = 1.0 - (ranks / len(ranks))
        
        # Non-linear boost: Make the top 20% really stand out
        # If stability > 0.8, score stays high. If < 0.5, drops fast.
        stability = stability ** 3 
        
        return stability

    def calculate_hk_scale(self, X, Y):
        """
        Calculates S using 'Robust Housekeeping' genes.
        Distinguishes Technical Gain (All genes shift) from Bio Bias (Only markers shift).
        """
        # 1. Determine Noise Floor (Avoid zero-variance trap)
        global_min = np.min(X)
        noise_buffer = 1e-5 if np.max(X) < 100 else 1.0
        dropout_threshold = global_min + noise_buffer
        
        # 2. Identify Expressed Genes
        mu_x = np.mean(X, axis=1)
        sigma_x = np.std(X, axis=1)
        is_expressed = mu_x > dropout_threshold
        
        if np.sum(is_expressed) < 10:
            return 1.0 # Fallback
            
        valid_indices = np.where(is_expressed)[0]
        
        # 3. Select Stable Anchors (Bottom 20% CoV)
        cov = sigma_x[valid_indices] / (np.abs(mu_x[valid_indices]) + EPS)
        
        # Hyperparam: hk_percentile (default 0.2)
        k_pct = getattr(self.hp, 'hk_percentile', 0.2)
        k = max(5, int(len(cov) * k_pct))
        
        best_sub_indices = np.argsort(cov)[:k]
        stable_idx = valid_indices[best_sub_indices]
        
        # 4. Calculate Ratio on Anchors Only
        X_hk = X[stable_idx, :]
        Y_hk = Y[stable_idx, :]
        
        mx = np.median(X_hk[X_hk > dropout_threshold])
        if np.isnan(mx): mx = 1.0
            
        my = np.median(Y_hk[Y_hk > dropout_threshold])
        if np.isnan(my) or my < EPS: my = 1.0
        
        S_hk = mx / (my + EPS)
        
        # Safety Clamps (Prevent explosion on bad data)
        S_hk = np.clip(S_hk, 0.01, 100.0)
        
        if self.debug:
            print(f"  Robust HK Anchor: Used {k} stable genes.")
            print(f"  Calculated S_HK = {S_hk:.4f}")
            
        return S_hk

    def adaptive_preprocessing(self, X, Y):
        """
        v5.0: Pure Local Preprocessing.
        We do NOT apply global scaling here because ComBat (the Prior) 
        handles the global scale estimation for us.
        """
        # Pre-Calculate Anchor Scores (for Local Trust)
        self.gene_stability_scores = self.calculate_gene_stability(X)
        
        # Transform without global scaling - ComBat priors handle scale
        X_prime = np.arcsinh(X)
        Y_prime = np.arcsinh(Y)
        
        return X_prime, Y_prime

    def pathway_execution(self, X_prime, Y_prime, pathway_indices, C_null, global_gene_idxs, alpha_prior_vec, beta_prior_vec, tau_override=None, top_k_override=None):
        """
        Execute peer finding and estimation for one pathway with gene-specific priors (v5.0).
        
        Args:
            alpha_prior_vec: Vector of gene-specific alpha priors (K,) for this pathway
            beta_prior_vec: Vector of gene-specific beta priors (K,) for this pathway
        """
        # Use override parameters if provided (for two-phase operation)
        tau = tau_override if tau_override is not None else self.hp.tau
        top_k_percent = top_k_override if top_k_override is not None else self.hp.top_k_percent
        
        # 1. Navigation (Centered Cosine on Raw)
        X_k = X_prime[pathway_indices, :]
        Y_k = Y_prime[pathway_indices, :]
        
        # Shapes match? (Correlation 1.0 means shapes are identical)
        K_raw = centered_cosine_similarity(X_k, Y_k)
        
        # Top-K Gating (with dynamic parameters)
        L_raw = tau * K_raw
        
        N_ref = L_raw.shape[0]
        k_neighbors = max(5, int(N_ref * top_k_percent))
        k_cutoff_idx = N_ref - k_neighbors
        
        # Hard Attention Mask (Top-K only)
        mask = np.full_like(L_raw, -np.inf)
        for col in range(L_raw.shape[1]):
            cutoff_val = np.partition(L_raw[:, col], k_cutoff_idx)[k_cutoff_idx]
            keep_mask = L_raw[:, col] >= cutoff_val
            mask[keep_mask, col] = L_raw[keep_mask, col]
        
        # Softmax & Weights
        null_row = np.full((1, K_raw.shape[1]), tau * C_null)
        L_aug = np.vstack([mask, null_row])
        P_aug = softmax(L_aug, axis=0)
        P_yx = P_aug[:-1, :]
        P_null = P_aug[-1, :]
        
        w_y = (1 - P_null)
        w_y_norm = w_y / (np.sum(w_y) + EPS)
        w_x = P_yx @ w_y_norm
        w_x_norm = w_x / (np.sum(w_x) + EPS)
        
        # 2. Local Estimation (Method of Moments on Raw Data)
        mu_x, var_x = weighted_stats(X_k, w_x_norm)
        mu_y, var_y = weighted_stats(Y_k, w_y_norm)
        
        # 3. Trust Gate (Pure Local Trust)
        X_virtual = X_k @ P_yx
        X_c = X_virtual - np.sum(X_virtual * w_y_norm, axis=1)[:, None]
        Y_c = Y_k - np.sum(Y_k * w_y_norm, axis=1)[:, None]
        
        cov = np.sum(w_y_norm * X_c * Y_c, axis=1)
        denom = np.sqrt(np.sum(w_y_norm*X_c**2, axis=1) * np.sum(w_y_norm*Y_c**2, axis=1))
        rho = cov / (denom + EPS)
        
        base_trust = np.maximum(0, rho)**2
        stability = self.gene_stability_scores[global_gene_idxs]
        
        # Logic: Trust Correlation OR Stability (Anchors)
        effective_trust = np.maximum(base_trust, stability)
        
        # 4. Adaptive Shrinkage (Gene-Specific Priors)
        alpha_raw = np.sqrt((var_x + EPS) / (var_y + EPS))
        beta_raw = mu_x - alpha_raw * mu_y
        
        # Interpolate with gene-specific priors:
        # Trust = 1.0 -> Use Local (Preserve Biology)
        # Trust = 0.0 -> Use Prior (Remove Artifact)
        alpha_est = alpha_prior_vec + effective_trust * (alpha_raw - alpha_prior_vec)
        beta_est = beta_prior_vec + effective_trust * (beta_raw - beta_prior_vec)
        
        # Fidelity metric
        H_norm = safe_entropy(P_yx) / np.log(N_ref + EPS)
        omega = np.mean((1 - H_norm) * (1 - P_null))
        
        # Metrics
        metrics = {
            'avg_trust': np.mean(effective_trust),
            'avg_alpha': np.mean(alpha_est),
            'avg_beta': np.mean(beta_est),
            'avg_correlation': np.mean(rho),
            'avg_stability': np.mean(stability)
        }
        
        return omega, alpha_est, beta_est, K_raw.flatten(), metrics
        k_cutoff_idx = N_ref - k_neighbors
        
        mask = np.full_like(L_raw, -np.inf)
        for col in range(L_raw.shape[1]):
            cutoff_val = np.partition(L_raw[:, col], k_cutoff_idx)[k_cutoff_idx]
            keep_mask = L_raw[:, col] >= cutoff_val
            mask[keep_mask, col] = L_raw[keep_mask, col]
            
        # Softmax & Weights
        null_row = np.full((1, K_raw.shape[1]), tau * C_null)
        L_aug = np.vstack([mask, null_row])
        P_aug = softmax(L_aug, axis=0)
        P_yx = P_aug[:-1, :]
        P_null = P_aug[-1, :]
        
        w_y = (1 - P_null)
        w_y_norm = w_y / (np.sum(w_y) + EPS)
        w_x = P_yx @ w_y_norm
        w_x_norm = w_x / (np.sum(w_x) + EPS)
        
        # 2. Estimation (Method of Moments)
        mu_x = np.sum(X_k * w_x_norm, axis=1)
        var_x = np.sum(w_x_norm * (X_k - mu_x[:, None])**2, axis=1)
        mu_y = np.sum(Y_k * w_y_norm, axis=1)
        var_y = np.sum(w_y_norm * (Y_k - mu_y[:, None])**2, axis=1)
        
        # 3. Clean Trust Logic (v4.2 - Max-Pooling)
        
        # A. Calculate Base Trust (Correlation of Projected Peers)
        X_virtual = X_k @ P_yx
        X_c = X_virtual - np.sum(X_virtual * w_y_norm, axis=1)[:, None]
        Y_c = Y_k - np.sum(Y_k * w_y_norm, axis=1)[:, None]
        cov = np.sum(w_y_norm * X_c * Y_c, axis=1)
        denom = np.sqrt(np.sum(w_y_norm*X_c**2, axis=1) * np.sum(w_y_norm*Y_c**2, axis=1))
        rho = cov / (denom + EPS)
        
        correlation_trust = np.maximum(0, rho)**2
        
        # B. Hybrid Trust Logic (Max-pooling)
        # Get stability scores for these specific genes
        local_stability = self.gene_stability_scores[global_gene_idxs]
        
        # Pure max-pooling trust calculation
        effective_trust = np.maximum(correlation_trust, local_stability)
        
        # 4. Raw Correction Calculation
        alpha_raw = (np.sqrt(var_x) + EPS) / (np.sqrt(var_y) + EPS)
        beta_raw = mu_x - alpha_raw * mu_y
        
        # 5. Adaptive Shrinkage toward Learned Prior (v4.2 Key Innovation)
        # If Trust=1, use Raw. If Trust=0, use Prior.
        alpha_prop = alpha_prior + effective_trust * (alpha_raw - alpha_prior)
        beta_prop = beta_prior + effective_trust * (beta_raw - beta_prior)
        
        # Metrics
        def gini(w):
            w_s = np.sort(w)
            return 1.0 - 2.0 * np.sum((np.cumsum(w_s) - w_s/2)/np.sum(w_s))/len(w)
            
        metrics = {
            'gini_x': gini(w_x_norm),
            'neff_x': 1.0/(np.sum(w_x_norm**2)+EPS),
            'avg_trust': np.mean(effective_trust),
            'avg_alpha': np.mean(alpha_prop),
            'alpha_raw': alpha_raw,
            'beta_raw': beta_raw,
            'alpha_prior': alpha_prior,
            'beta_prior': beta_prior
        }
        
        # Omega (Pathway Fidelity)
        H_norm = safe_entropy(P_yx) / np.log(N_ref + EPS)
        omega = np.mean((1 - H_norm) * (1 - P_null))
        
        return omega, alpha_prop, beta_prop, K_raw.flatten(), metrics

    def align(self, ref_data: BatchData, target_data: BatchData) -> Tuple[np.ndarray, Dict]:
        """
        Execute POSSE v5.0 alignment with ComBat-Initialized gene-specific priors.
        """
        if self.debug: print("Starting POSSE v5.0 (ComBat-Initialized)...")
        
        # Intersection
        common, idx_x, idx_y = np.intersect1d(
            ref_data.gene_indices, target_data.gene_indices, return_indices=True
        )
        X = ref_data.data[idx_x]
        Y = target_data.data[idx_y]
        
        # Preprocessing (Pure Local)
        X_p, Y_p = self.adaptive_preprocessing(X, Y)
        if self.debug: print(f"  Pure local preprocessing applied")
        
        # 1. INITIALIZE PRIORS WITH COMBAT
        # This gives us the best "Blind" guess for artifacts
        if self.debug: print("  Calculating Initial ComBat Priors...")
        alpha_prior, beta_prior = self.combat.compute_baseline(X_p, Y_p)
        
        if self.debug:
            print(f"  ComBat priors: Alpha mean={np.mean(alpha_prior):.3f}, Beta mean={np.mean(beta_prior):.3f}")
        
        # Prepare Pathways
        gene_map = {g: i for i, g in enumerate(common)}
        pathway_idxs = []
        for name, genes in self.pathway_dict.items():
            idxs = [gene_map[g] for g in genes if g in gene_map]
            if len(idxs) >= self.hp.min_pathway_size:
                pathway_idxs.append(np.array(idxs))
                
        if not pathway_idxs:
            raise RuntimeError("CRITICAL: No valid pathways found for POSSE correction. "
                             f"Loaded {len(self.pathway_dict)} pathways but none have sufficient genes "
                             f"(min_pathway_size={self.hp.min_pathway_size}) that match your data. "
                             f"Check gene symbol compatibility or pathway database.")
                
        G = len(common)
        C_null = 0.0
        
        # Diagnostic Collectors
        diag_history = []
        
        # Iterative Refinement
        for t in range(self.hp.max_iter):
            if self.debug: 
                print(f"  Iter {t+1}: Refining Priors...")
                print(f"    Current Alpha range: [{np.min(alpha_prior):.3f}, {np.max(alpha_prior):.3f}]")
                print(f"    Current Beta range: [{np.min(beta_prior):.3f}, {np.max(beta_prior):.3f}]")
            
            sig_a = np.zeros(G)
            sig_b = np.zeros(G)
            sig_w = np.zeros(G) + EPS  # Prevent div/0
            
            K_vals = []
            metrics_acc = {'trust': [], 'alpha': [], 'beta': [], 'correlation': [], 'stability': []}
            
            for p_idx in pathway_idxs:
                # Pass the GENE-SPECIFIC priors for this pathway
                omega, a_p, b_p, k_raw, mets = self.pathway_execution(
                    X_p, Y_p, p_idx, C_null, p_idx,
                    alpha_prior[p_idx], beta_prior[p_idx]  # Vector inputs
                )
                
                K_vals.extend(k_raw)
                metrics_acc['trust'].append(mets['avg_trust'])
                metrics_acc['alpha'].append(mets['avg_alpha'])
                metrics_acc['beta'].append(mets['avg_beta'])
                metrics_acc['correlation'].append(mets['avg_correlation'])
                metrics_acc['stability'].append(mets['avg_stability'])
                
                # Accumulate Votes
                # We weight votes by Pathway Fidelity (Omega)
                for i, g_idx in enumerate(p_idx):
                    sig_a[g_idx] += a_p[i] * omega
                    sig_b[g_idx] += b_p[i] * omega
                    sig_w[g_idx] += omega
            
            # UPDATE PRIORS (Consensus)
            # The new prior is the consensus of the local experts
            # Genes with no pathway coverage fall back to their previous prior (ComBat)
            mask_covered = sig_w > EPS
            
            new_alpha = alpha_prior.copy()
            new_beta = beta_prior.copy()
            
            new_alpha[mask_covered] = sig_a[mask_covered] / sig_w[mask_covered]
            new_beta[mask_covered] = sig_b[mask_covered] / sig_w[mask_covered]
            
            # Smooth update to prevent oscillation
            alpha_prior = alpha_prior * 0.2 + new_alpha * 0.8
            beta_prior = beta_prior * 0.2 + new_beta * 0.8
            
            # Update Null Threshold
            if K_vals:
                C_null = (1 - self.hp.eta)*C_null + self.hp.eta*np.percentile(K_vals, 5)
            
            # Diagnostic Snapshot
            diag_history.append({
                'iter': t,
                'avg_trust': np.mean(metrics_acc['trust']),
                'avg_alpha': np.mean(metrics_acc['alpha']),
                'avg_beta': np.mean(metrics_acc['beta']),
                'avg_correlation': np.mean(metrics_acc['correlation']),
                'avg_stability': np.mean(metrics_acc['stability']),
                'c_null': C_null,
                'genes_covered': np.sum(mask_covered),
                'alpha_range': [np.min(alpha_prior), np.max(alpha_prior)],
                'beta_range': [np.min(beta_prior), np.max(beta_prior)]
            })
            
            if self.debug:
                print(f"    Genes covered by pathways: {np.sum(mask_covered)}/{G}")
                print(f"    Avg trust: {np.mean(metrics_acc['trust']):.3f}")
                print(f"    Avg correlation: {np.mean(metrics_acc['correlation']):.3f}")
        
        # Final Application
        Y_corr = alpha_prior[:, None] * Y_p + beta_prior[:, None]
        Y_final = np.sinh(Y_corr)
        
        # Metadata
        metadata = {
            'alpha_final': alpha_prior,
            'beta_final': beta_prior,
            'alpha_final_mean': np.mean(alpha_prior),
            'beta_final_mean': np.mean(beta_prior),
            'S_diagnostic': 1.0,  # v5.0 doesn't use global scaling
            'genes_covered': np.sum(mask_covered),
            'genes_total': G,
            'diagnostics': diag_history,
            'version': '5.0'
        }
        
        if self.debug:
            print(f"POSSE v5.0 correction completed. Output shape: {Y_final.shape}")
            print(f"Final Alpha mean: {metadata['alpha_final_mean']:.4f}")
            print(f"Final Beta mean: {metadata['beta_final_mean']:.4f}")
            print(f"Genes covered: {metadata['genes_covered']}/{metadata['genes_total']}")
        
        # Reassemble logic
        return self._reassemble(target_data, common, Y_final, alpha_prior, beta_prior, metadata)
    
    def _reassemble(self, target_data, common_genes, corrected_common, alphas, betas, metadata):
        """Reassemble corrected data back to original target format"""
        Y_final = target_data.data.copy()
        t_map = {g: i for i, g in enumerate(target_data.gene_indices)}
        
        for i, gene in enumerate(common_genes):
            if gene in t_map:
                Y_final[t_map[gene]] = corrected_common[i]
        
        # Add per-gene parameters to metadata for diagnostic output
        metadata['alpha_final'] = alphas
        metadata['beta_final'] = betas
        metadata['common_genes'] = common_genes
        
        return Y_final, metadata
        
        a_final = np.ones(G)  # Default to identity
        b_final = np.zeros(G)
        
        # Apply POSSE corrections where pathways provided coverage
        covered_mask = ~uncovered_mask
        if np.any(covered_mask):
            a_final[covered_mask] = sig_a[covered_mask] / (sig_w[covered_mask] + EPS)
            b_final[covered_mask] = sig_b[covered_mask] / (sig_w[covered_mask] + EPS)
        
        if self.debug and np.any(uncovered_mask):
            n_uncovered = np.sum(uncovered_mask)
            print(f"  {n_uncovered}/{G} genes had no pathway coverage (using identity transform)")
        
        # Apply to common genes
        Y_corr = a_final[:, None] * Y_p + b_final[:, None]
        Y_common_final = np.sinh(Y_corr)
        
        # Reassemble full target data
        Y_final = self._reassemble(target_data, common, Y_common_final, a_final, b_final)
        
        return Y_final, {
            "S_diagnostic": S_diag, 
            "diagnostics": diag_history,
            "alpha_final_mean": np.mean(a_final),
            "beta_final_mean": np.mean(b_final),
            "genes_covered": np.sum(covered_mask),
            "genes_total": G
        }
