import numpy as np
from scipy.special import softmax
from scipy.stats import pearsonr
from dataclasses import dataclass
from typing import List, Dict, Tuple, Optional, Union
import os
import sys

try:
    import gseapy as gp
except ImportError:
    gp = None

try:
    import pandas as pd
except ImportError:
    pd = None

# ==========================================
# Global Constants (Derived from System)
# ==========================================
# Replaces hardcoded 1e-9 with machine precision
EPS = np.finfo(float).eps 

# ==========================================
# Utility Functions
# ==========================================

def weighted_pearson_corr(x: np.ndarray, y: np.ndarray, w: np.ndarray) -> float:
    """
    Compute weighted Pearson correlation between two vectors of equal length.
    w: Weights for the observations (samples).
    """
    if len(x) != len(y):
        return 0.0
        
    w_sum = np.sum(w)
    if w_sum == 0:
        return 0.0
        
    # Weighted Means
    mu_x = np.sum(x * w) / w_sum
    mu_y = np.sum(y * w) / w_sum
    
    # Centered vectors
    x_c = x - mu_x
    y_c = y - mu_y
    
    # Weighted Covariance
    cov = np.sum(w * x_c * y_c) / w_sum
    
    # Weighted Variances
    var_x = np.sum(w * x_c**2) / w_sum
    var_y = np.sum(w * y_c**2) / w_sum
    
    if var_x < EPS or var_y < EPS:
        return 0.0
        
    return cov / np.sqrt(var_x * var_y)

def safe_entropy(p: np.ndarray) -> np.ndarray:
    """Compute entropy safely handling p=0 cases."""
    # p shape: (N_ref, M_target)
    # We want entropy per Target sample (columns)
    return -np.sum(p * np.log(p + EPS), axis=0)

def cosine_similarity(X: np.ndarray, Y: np.ndarray) -> np.ndarray:
    """
    Compute cosine similarity between reference (X) and target (Y) samples.
    
    Args:
        X: Reference data (genes x samples_ref)
        Y: Target data (genes x samples_target)
    
    Returns:
        Similarity matrix (samples_ref x samples_target)
    """
    # Normalize vectors to unit length
    X_norm = X / (np.linalg.norm(X, axis=0, keepdims=True) + EPS)
    Y_norm = Y / (np.linalg.norm(Y, axis=0, keepdims=True) + EPS)
    
    # Compute cosine similarity: X^T @ Y
    return X_norm.T @ Y_norm

def centered_cosine_similarity(X: np.ndarray, Y: np.ndarray) -> np.ndarray:
    """
    Compute centered cosine similarity (Pearson correlation) between samples.
    X: (Genes, N)
    Y: (Genes, M)
    Returns: (N, M) matrix where (i,j) is sim(Ref_i, Target_j)
    """
    # Center the data (subtract mean across genes for each sample)
    X_mean = np.mean(X, axis=0, keepdims=True)
    Y_mean = np.mean(Y, axis=0, keepdims=True)
    
    X_centered = X - X_mean
    Y_centered = Y - Y_mean
    
    # Compute norms
    X_norms = np.linalg.norm(X_centered, axis=0, keepdims=True) + EPS
    Y_norms = np.linalg.norm(Y_centered, axis=0, keepdims=True) + EPS
    
    # Normalize and compute similarity matrix
    X_norm = X_centered / X_norms
    Y_norm = Y_centered / Y_norms
    
    # Return (N, M) similarity matrix
    return X_norm.T @ Y_norm

def rbf_kernel_similarity(X: np.ndarray, Y: np.ndarray, gamma: float = 0.1) -> np.ndarray:
    """
    Compute RBF Kernel similarity (Gaussian): exp(-gamma * ||x - y||^2)
    
    Args:
        X: Reference data (genes x samples_ref)
        Y: Target data (genes x samples_target)  
        gamma: RBF kernel parameter (higher = stricter matching)
    
    Returns:
        Similarity matrix (samples_ref x samples_target)
    """
    # Efficient Euclidean Distance Matrix Calculation:
    # ||x - y||^2 = ||x||^2 + ||y||^2 - 2<x, y>
    
    X_sq = np.sum(X**2, axis=0)  # Shape (N,)
    Y_sq = np.sum(Y**2, axis=0)  # Shape (M,)
    XY = X.T @ Y                 # Shape (N, M)
    
    # Broadcast addition: (N,1) + (1,M) - (N,M)
    dists_sq = X_sq[:, np.newaxis] + Y_sq - 2*XY
    
    # Clip negative zeros (numerical precision issues)
    dists_sq = np.maximum(dists_sq, 0.0)
    
    # RBF Kernel: exp(-gamma * distance^2)
    return np.exp(-gamma * dists_sq)

# ==========================================
# Data Structures
# ==========================================

@dataclass
class BatchData:
    data: np.ndarray  # Shape (Genes, Samples)
    gene_indices: np.ndarray

@dataclass 
class PACEHyperparameters:
    tau: float = 10.0  # Temperature for softmax sharpness
    w_prior: float = 1.0  # Prior weight strength
    lambda_damp: float = 2.0  # Amplification dampener
    max_iter: int = 5  # Maximum refinement loops
    eta: float = 0.2  # Null threshold learning rate
    use_hard_gating: bool = True  # Enable Top-K masking to fix Simpson's Paradox
    hard_gating_ratio: float = 0.10  # Fraction of reference samples to keep (10% for aggressive minority matching)
    similarity_metric: str = 'centered_cosine'  # 'centered_cosine' or 'cosine' - similarity metric for reference matching
    use_v23_correction: bool = False  # Enable PACE v2.3 "correction on original" improvement
    use_v24_navigation: bool = False  # Enable PACE v2.4 "scale-decoupled navigation" improvement
    use_v25_rbf: bool = False  # Enable PACE v2.5 "RBF kernel navigation" improvement
    use_v30_activity: bool = False  # Enable PACE v3.0 "activity-gated consensus" improvement
    use_v31_pure_local: bool = False  # Enable PACE v3.1 "pure local estimation" improvement
    use_v32_iterative_consensus: bool = False  # Enable PACE v3.2 "iterative consensus" improvement
    use_v34_housekeeping_anchors: bool = False  # Enable PACE v3.4 "housekeeping anchors" improvement
    use_v35_affine_anchors: bool = False  # Enable PACE v3.5 "global affine anchors" improvement
    gamma: float = 0.5  # RBF kernel gamma parameter (for v2.5) / Activity sharpness (for v3.0)
    consensus_threshold: float = 0.3  # Variance threshold for consensus decision (for v3.2)
    hk_percentile: float = 0.20  # Use bottom 20% variance genes as anchors (for v3.4/v3.5)

# ==========================================
# ComBat Baseline Implementation
# ==========================================

class ComBatBaseline:
    """Standard ComBat implementation for global baseline."""
    
    def compute_baseline(self, X_prime: np.ndarray, Y_prime: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """Compute global ComBat baseline parameters using Empirical Bayes."""
        G, N = X_prime.shape
        _, M = Y_prime.shape
        
        # Compute sample means and variances
        mu_X = np.mean(X_prime, axis=1)
        mu_Y = np.mean(Y_prime, axis=1)
        
        var_X = np.var(X_prime, axis=1, ddof=1)
        var_Y = np.var(Y_prime, axis=1, ddof=1)
        
        # Raw parameter estimates
        alpha_raw = np.sqrt((var_X + EPS) / (var_Y + EPS))
        beta_raw = mu_X - alpha_raw * mu_Y
        
        # Empirical Bayes shrinkage for alpha (scale parameters)
        # Inverse Gamma prior for variance
        V_g = alpha_raw**2
        mean_V = np.mean(V_g)
        var_V = np.var(V_g) + EPS
        
        # Method of moments for Inverse Gamma parameters
        lambda_prior = (mean_V**2 / var_V) + 2
        theta_prior = mean_V * (lambda_prior - 1)
        
        # Posterior parameters
        alpha_sq_post = (theta_prior + 0.5 * M * V_g) / (lambda_prior + 0.5 * M)
        alpha_glob = np.sqrt(alpha_sq_post)
        
        # Normal prior for beta (location parameters)
        gamma_prior = np.mean(beta_raw)
        tau_sq_prior = np.var(beta_raw) + EPS
        
        # Posterior mean
        precision_prior = 1.0 / tau_sq_prior
        precision_data = M / var_Y
        
        beta_glob = (precision_prior * gamma_prior + precision_data * beta_raw) / (precision_prior + precision_data)
        
        return alpha_glob, beta_glob 

# ==========================================
# Pathway Data Functions
# ==========================================

def read_gmt(filename: str) -> Dict[str, List[str]]:
    pathways = {}
    if not os.path.exists(filename):
        return pathways
    try:
        with open(filename, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 3:
                    term = parts[0]
                    genes = [g for g in parts[2:] if g]
                    pathways[term] = genes
    except Exception as e:
        print(f"Error reading GMT: {e}")
    return pathways

def load_pathways(name: str, organism: str = 'Human', save: Optional[str] = None) -> Dict[str, List[str]]:
    """Simplified loader with sensible defaults."""
    target_file = name if (name.endswith('.gmt') or name.endswith('.txt')) else save
    
    if target_file and os.path.exists(target_file):
        return read_gmt(target_file)
        
    try:
        print(f"Downloading {name}...")
        pathways = gp.get_library(name=name, organism=organism)
        if save:
            pass # Save logic omitted for brevity
        return pathways
    except Exception:
        return {}

# ==========================================
# PACE v2.2 Implementation
# ==========================================

class PACE_v22:
    """
    Pathway-Adaptive Consensus Estimator (Robust Hybrid) v2.2
    
    Implements the modular PACE algorithm with:
    1. Adaptive pre-processing with global scale factor
    2. Iterative pathway-based parameter estimation
    3. Dynamic null threshold learning
    4. Orphan gene handling with lookup table extrapolation
    """
    
    def __init__(self, 
                 pathway_dict: Dict[str, List[str]] = None,
                 pathway_source: str = 'hallmark',
                 organism: str = 'Human',
                 save_path: str = None,
                 hyperparams: PACEHyperparameters = None):
        
        self.hyperparams = hyperparams or PACEHyperparameters()
        
        # Pathway Loading Logic
        if pathway_dict:
            self.pathway_dict = pathway_dict
        else:
            name = 'MSigDB_Hallmark_2020' if pathway_source == 'hallmark' else pathway_source
            self.pathway_dict = load_pathways(name, organism, save=save_path)
        
        self.combat_baseline = ComBatBaseline()
        
        # Initialize coherence scores (will be calculated during alignment)
        self.coherence_scores = None
    
    def adaptive_preprocessing(self, X: np.ndarray, Y: np.ndarray) -> Tuple[np.ndarray, np.ndarray, float]:
        """
        Step 2: Adaptive Pre-Processing
        Compute global scale factor and apply VSN-like transformation.
        
        For PACE v3.0, returns additional Y_nav (unscaled navigation data)
        For PACE v3.1, disables global scaling entirely (S=1.0)
        For PACE v3.2, uses iterative consensus to discover global artifacts
        For PACE v3.4, uses housekeeping anchors to distinguish biological vs technical
        For PACE v3.5, uses housekeeping affine anchors for shift + scale correction
        """
        if self.hyperparams.use_v35_affine_anchors:
            # PACE v3.5: Global Affine Anchors - Use stable genes for shift + scale
            alpha_glob, beta_glob = self.calculate_hk_affine(X, Y)
            
            # Apply Affine Transform: Y_corrected = (Y - beta) * alpha
            Y_shifted = Y - beta_glob  # Remove additive bias first
            Y_scaled = Y_shifted * alpha_glob  # Apply multiplicative correction
            
            # Apply arcsinh for numerical stability
            X_prime = np.arcsinh(X)
            Y_prime = np.arcsinh(Y_scaled)
            Y_raw = np.arcsinh(Y)  # Keep unscaled for navigation safety
            
            print(f"🔗 PACE v3.5: Global Affine Anchors - α={alpha_glob:.4f}, β={beta_glob:.4f}")
            return X_prime, Y_prime, Y_raw, (alpha_glob, beta_glob)
        elif self.hyperparams.use_v34_housekeeping_anchors:
            # PACE v3.4: Housekeeping Anchors - Use stable genes to calculate S
            S = self.calculate_hk_scale(X, Y)
            X_prime = np.arcsinh(X)
            Y_prime = np.arcsinh(Y * S)  # Apply HK correction globally
            Y_raw = np.arcsinh(Y)  # Keep unscaled for navigation safety
            print(f"🔗 PACE v3.4: Housekeeping Anchors - S_HK={S:.4f}")
            return X_prime, Y_prime, Y_raw, S
        elif self.hyperparams.use_v32_iterative_consensus:
            # PACE v3.2: Iterative Consensus - Start with S=1.0, let pathways vote
            X_prime = np.arcsinh(X)
            Y_prime = np.arcsinh(Y)
            S = 1.0  # Initial S, will be updated by consensus
            print(f"🔬 PACE v3.2: Iterative Consensus - Starting with S=1.0, pathways will vote")
            return X_prime, Y_prime, S
        elif self.hyperparams.use_v31_pure_local:
            # PACE v3.1: Pure Local Estimation - NO Global Scaling (S=1.0)
            # We perform arcsinh transform only for numerical stability
            X_prime = np.arcsinh(X)
            Y_prime = np.arcsinh(Y)
            S = 1.0  # No global scaling
            print(f"🔬 PACE v3.1: Pure Local Estimation - Global scaling disabled (S=1.0)")
            return X_prime, Y_prime, S
        
        # Standard global scale factor calculation for other versions
        # Global Scale Factor (S)
        x_positive = X[X > 0]
        y_positive = Y[Y > 0]
        
        median_x = np.median(x_positive) if len(x_positive) > 0 else 1.0
        median_y = np.median(y_positive) if len(y_positive) > 0 else 1.0
        
        # Robust scaling factor calculation
        EPS = 1e-8
        if median_x < EPS: 
            median_x = 1.0
        if median_y < EPS: 
            median_y = 1.0
            
        S = median_x / median_y
        
        # Absolute safety clamp: Prevent S from being extreme
        S = np.clip(S, 0.1, 10.0)
        
        # VSN-like transformation
        X_prime = np.arcsinh(X)
        
        if self.hyperparams.use_v30_activity:
            # PACE v3.0: Return both scaled and unscaled target data
            Y_nav = np.arcsinh(Y)  # Unscaled for shape navigation
            Y_prime = np.arcsinh(Y * S)  # Scaled for activity calculation & correction
            return X_prime, Y_prime, Y_nav, S
        else:
            # Standard PACE: Only scaled data
            Y_prime = np.arcsinh(Y * S)  # Scaled for correction
            return X_prime, Y_prime, S
    
    def calculate_global_reference_mask(self, X_prime: np.ndarray, Y_prime: np.ndarray) -> np.ndarray:
        """
        Calculate global Top-K reference mask to use across all pathways.
        This fixes the Simpson's Paradox issue by ensuring only a small subset of reference samples
        are used globally, rather than different subsets per pathway.
        """
        N, M = X_prime.shape[1], Y_prime.shape[1]
        
        if not self.hyperparams.use_hard_gating:
            return None
            
        k_global = max(10, int(N * self.hyperparams.hard_gating_ratio))
        k_global = min(k_global, N)
        
        # Calculate global similarity across all genes using selected metric
        if self.hyperparams.similarity_metric == 'centered_cosine':
            K_global = centered_cosine_similarity(X_prime, Y_prime)  # (N, M)
        elif self.hyperparams.similarity_metric == 'cosine':
            K_global = cosine_similarity(X_prime, Y_prime)  # (N, M)
        else:
            # Default to centered cosine
            K_global = centered_cosine_similarity(X_prime, Y_prime)  # (N, M)
        
        # For each target sample, find the top-K reference samples globally
        global_ref_mask = np.zeros((N, M), dtype=bool)
        for m in range(M):
            # Get top k_global reference samples for this target
            top_k_indices = np.argpartition(K_global[:, m], -k_global)[-k_global:]
            global_ref_mask[top_k_indices, m] = True
        
        avg_refs_per_target = np.mean(np.sum(global_ref_mask, axis=0))
        return global_ref_mask
    
    def calculate_hk_scale(self, X: np.ndarray, Y: np.ndarray) -> float:
        """
        Calculates S using 'Robust Housekeeping' genes.
        
        Fixes the 'Zero-Variance Trap':
        - Genes that are all zeros (or near dataset minimum) are BANNED.
        - Only genes clearly distinguishable from the noise floor are eligible.
        """
        # 1. Determine the Noise Floor
        # We look at the entire reference matrix to find the true bottom
        global_min = np.min(X)
        
        # Define a safe deviation away from the minimum
        # If min is 0, this requires genes to have mean > 0.01
        # If min is negative (centered), it requires them to be above the floor
        noise_buffer = 1e-5 
        if np.max(np.abs(X)) > 100: # Heuristic: If counts, buffer is 1.0
            noise_buffer = 1.0
        
        dropout_threshold = global_min + noise_buffer
        
        # 2. Filter Eligible Genes (Must be expressed)
        # We use Mean to be robust against outliers
        mu_x = np.mean(X, axis=1)
        sigma_x = np.std(X, axis=1)
        
        # THE GATE: Gene must be "On" to be stable
        # This prevents selecting genes that are 0.0 everywhere
        is_expressed = mu_x > dropout_threshold
        
        if np.sum(is_expressed) < 10:
            print("  [WARNING] Almost no genes expressed above noise floor! Falling back to global median.")
            valid_indices = np.arange(len(mu_x)) # Fallback to all
        else:
            valid_indices = np.where(is_expressed)[0]
            
        # 3. Select Stable Genes (Low CoV) from the Expressed set
        # CoV = Sigma / Mu
        # Add epsilon to Mu to prevent div-by-zero for edge cases
        cov = sigma_x[valid_indices] / (mu_x[valid_indices] + EPS)
        
        # Select bottom K percentile (Most stable)
        k = int(len(cov) * self.hyperparams.hk_percentile)
        k = max(5, k) # Ensure at least 5 genes
        
        # Sort by CoV ascending (Low variance first)
        best_sub_indices = np.argsort(cov)[:k]
        stable_idx = valid_indices[best_sub_indices]
        
        # 4. Calculate S on Stable Genes ONLY
        X_hk = X[stable_idx, :]
        Y_hk = Y[stable_idx, :]
        
        # Calculate scaling factor
        # We check Y for zeros too, just in case
        mx = np.median(X_hk[X_hk > dropout_threshold]) 
        if np.isnan(mx): mx = 1.0
            
        my = np.median(Y_hk[Y_hk > dropout_threshold])
        if np.isnan(my) or my < EPS: my = 1.0
        
        S_hk = np.clip(mx / (my + EPS), 0.1, 10.0)
        
        print(f"  🔗 Robust HK Anchor: Found {len(valid_indices)} expressed genes.")
        print(f"  🔗 Selected {k} stable anchors (Min Mean: {np.min(mu_x[stable_idx]):.4f}).")
        print(f"  🔗 Calculated S_HK = {S_hk:.4f}")
        
        return S_hk
    
    def calculate_hk_affine(self, X: np.ndarray, Y: np.ndarray) -> Tuple[float, float]:
        """
        Calculates Global Shift (Beta) and Scale (Alpha) using HK genes.
        Robustly aligns the 'inert' biological structure using affine transformation.
        
        Returns:
            Tuple of (alpha_glob, beta_glob) for affine transform: Y_corrected = (Y - beta) * alpha
        """
        # 1. Identify Stable Genes in Reference (Low CoV)
        # Using Median/IQR for robustness against outliers
        mu_x = np.median(X, axis=1)
        # IQR as robust sigma
        sigma_x = np.percentile(X, 75, axis=1) - np.percentile(X, 25, axis=1)
        cov = sigma_x / (np.abs(mu_x) + EPS)
        
        # Filter for expressed genes
        global_min = np.min(X)
        mask_expr = mu_x > (global_min + 0.1)
        if np.sum(mask_expr) < 10:
            valid_idx = np.arange(len(mu_x))
        else:
            valid_idx = np.where(mask_expr)[0]
            
        # Select bottom K% CoV
        cov = cov[valid_idx]
        k = max(5, int(len(cov) * self.hyperparams.hk_percentile))
        best_idx = np.argsort(cov)[:k]
        stable_idx = valid_idx[best_idx]
        
        # 2. Calculate Affine Parameters on HK Genes
        X_hk = X[stable_idx, :]
        Y_hk = Y[stable_idx, :]
        
        # Shift (Median Difference)
        med_x = np.median(X_hk)
        med_y = np.median(Y_hk)
        beta_glob = med_y - med_x
        
        # Scale (IQR Ratio)
        iqr_x = np.percentile(X_hk, 75) - np.percentile(X_hk, 25)
        iqr_y = np.percentile(Y_hk, 75) - np.percentile(Y_hk, 25)
        alpha_glob = iqr_x / (iqr_y + EPS)
        
        # Safety clamps
        alpha_glob = np.clip(alpha_glob, 0.5, 2.0)
        beta_glob = np.clip(beta_glob, -5.0, 5.0)
        
        print(f"  🔗 Robust HK Affine: Found {len(valid_idx)} expressed genes.")
        print(f"  🔗 Selected {k} stable anchors for affine calculation.")
        print(f"  🔗 Calculated α={alpha_glob:.4f}, β={beta_glob:.4f}")
        
        return alpha_glob, beta_glob
        print(f"  🔗 HK Affine: {len(stable_idx)} genes. Shift={beta_glob:.4f}, Scale={alpha_glob:.4f}")
        
        return alpha_glob, beta_glob
        """
        Calculate global Top-K reference mask to use across all pathways.
        This fixes the Simpson's Paradox issue by ensuring only a small subset of reference samples
        are used globally, rather than different subsets per pathway.
        """
        N, M = X_prime.shape[1], Y_prime.shape[1]
        
        if not self.hyperparams.use_hard_gating:
            return None
            
        k_global = max(10, int(N * self.hyperparams.hard_gating_ratio))
        k_global = min(k_global, N)
        
        print(f"  GLOBAL Top-K masking: selecting {k_global}/{N} ({100*k_global/N:.1f}%) reference samples for ALL pathways")
        
        # Calculate global similarity across all genes using selected metric
        if self.hyperparams.similarity_metric == 'centered_cosine':
            K_global = centered_cosine_similarity(X_prime, Y_prime)  # (N, M)
        elif self.hyperparams.similarity_metric == 'cosine':
            K_global = cosine_similarity(X_prime, Y_prime)  # (N, M)
        else:
            # Default to centered cosine
            K_global = centered_cosine_similarity(X_prime, Y_prime)  # (N, M)
        
        # For each target sample, find the top-K reference samples globally
        global_ref_mask = np.zeros((N, M), dtype=bool)
        for m in range(M):
            # Get top k_global reference samples for this target
            top_k_indices = np.argpartition(K_global[:, m], -k_global)[-k_global:]
            global_ref_mask[top_k_indices, m] = True
        
        avg_refs_per_target = np.mean(np.sum(global_ref_mask, axis=0))
        print(f"    Global masking: average {avg_refs_per_target:.1f} reference samples per target")
        
        return global_ref_mask

    def pathway_execution(self, X_prime: np.ndarray, Y_prime: np.ndarray, 
                         pathway_indices: np.ndarray, C_null: float, 
                         pathway_name: str = None, use_topk_masking: bool = True,
                         global_ref_mask: np.ndarray = None, Y_raw: np.ndarray = None) -> Tuple[float, np.ndarray, np.ndarray, np.ndarray, Dict[str, float]]:
        """
        Execute pathway-based parameter estimation with coherence weighting and Top-K masking.
        
        Args:
            Y_raw: Unscaled target data for v2.4 navigation (optional)
        
        Returns: fidelity, alphas, betas, raw_sims, metrics_dict
        """
        # Extract pathway-specific data
        X_k = X_prime[pathway_indices, :]  # (k_genes, N)
        Y_k = Y_prime[pathway_indices, :]  # (k_genes, M) - Scaled for correction
        
        # PACE v2.4: Use unscaled data for navigation if available
        if self.hyperparams.use_v24_navigation and Y_raw is not None:
            Y_k_nav = Y_raw[pathway_indices, :]  # (k_genes, M) - Unscaled for navigation
            print(f"      Using PACE v2.4: Scale-decoupled navigation")
        else:
            Y_k_nav = Y_k  # Use scaled data for both navigation and correction (v2.2/v2.3 behavior)
        
        # Local Navigation - Similarity calculation based on selected metric
        if self.hyperparams.use_v30_activity:
            # PACE v3.0: Activity-Gated Consensus (Shape + Activity)
            # PHASE 1: SHAPE NAVIGATION (Cosine)
            K_shape = centered_cosine_similarity(X_k, Y_k_nav)
            
            # PHASE 2: ACTIVITY GATING (1D Euclidean)
            # Calculate Pathway Activity Score (Scalar Mean)
            # Use Y_prime (Scaled) to match absolute levels against Reference
            act_x = np.mean(X_k, axis=0)  # Shape (N,)
            act_y = np.mean(Y_k, axis=0)  # Shape (M,) - Use scaled Y for activity
            
            # Calculate 1D Distance Matrix: |Ax - Ay|
            dist_act = np.abs(act_x[:, np.newaxis] - act_y[np.newaxis, :])
            
            # Convert to Similarity Gate (Gaussian)
            # Scale Gamma by expected range of activity (approx 1.0 for arcsinh)
            gamma = self.hyperparams.gamma
            K_act = np.exp(-gamma * (dist_act**2))
            
            # FUSE SIGNALS: Combined Similarity = Shape * Activity
            # If shapes match but levels don't, score -> 0.
            K_raw = K_shape * K_act
            
        elif self.hyperparams.use_v25_rbf:
            # PACE v2.5: RBF Kernel Navigation (Euclidean distance-based)
            n_genes = len(pathway_indices)
            scaled_gamma = self.hyperparams.gamma / (n_genes + EPS)
            K_raw = rbf_kernel_similarity(X_k, Y_k_nav, gamma=scaled_gamma)
            
        elif self.hyperparams.similarity_metric == 'centered_cosine':
            K_raw = centered_cosine_similarity(X_k, Y_k_nav)
        elif self.hyperparams.similarity_metric == 'cosine':
            K_raw = cosine_similarity(X_k, Y_k_nav)
        else:
            # Default to centered cosine
            K_raw = centered_cosine_similarity(X_k, Y_k_nav)
        
        # Calculate Similarity Logits
        if self.hyperparams.use_v25_rbf:
            # RBF kernel outputs similarities (0-1), not logits
            # We can still use tau for temperature scaling if desired
            tau = self.hyperparams.tau
            L_raw = tau * K_raw  # Optional temperature scaling
        else:
            # Standard cosine similarities converted to logits
            tau = self.hyperparams.tau
            L_raw = tau * K_raw
        
        # NEW: TOP-K MASKING LOGIC to fix Simpson's Paradox
        if use_topk_masking and global_ref_mask is None:
            # Per-pathway Top-K masking (original approach)
            N_ref = L_raw.shape[0]
            k_neighbors = max(5, int(N_ref * self.hyperparams.hard_gating_ratio))
            k_neighbors = min(k_neighbors, 100)  # Cap for efficiency
            
            # Generate mask (original logits for keep, -inf for mask)
            mask = np.full_like(L_raw, -np.inf)
            
            # This loop is unavoidable but fast for typical batch sizes
            for col in range(L_raw.shape[1]):
                # Get the top k_neighbors indices (handles ties correctly)
                top_k_indices = np.argpartition(L_raw[:, col], -k_neighbors)[-k_neighbors:]
                # Copy original logits for top-k samples only
                mask[top_k_indices, col] = L_raw[top_k_indices, col]
            
            L_masked = mask
        elif global_ref_mask is not None:
            # Global Top-K masking (new approach to fix Simpson's Paradox)
            # Apply global mask: set non-selected reference samples to -inf
            mask = np.full_like(L_raw, -np.inf)
            for col in range(L_raw.shape[1]):
                # Use only globally selected reference samples for this target
                selected_refs = global_ref_mask[:, col]
                mask[selected_refs, col] = L_raw[selected_refs, col]
            
            L_masked = mask
        else:
            L_masked = L_raw
        
        # Null-Augmented Softmax
        if self.hyperparams.use_v25_rbf:
            # RBF kernel special handling: similarities are already probabilities (0-1)
            # Null handling: If max similarity is low, trust drops
            max_sim = np.max(K_raw, axis=0)
            P_null = 1.0 - max_sim  # If best match is 0.9, Null is 0.1. If 0.1, Null is 0.9.
            
            # Normalize the masked similarities to sum to (1 - P_null)
            if use_topk_masking and global_ref_mask is None:
                # Use the masked similarities
                P_yx = L_masked / (tau + EPS)  # Convert back from logits if tau was applied
                P_yx = np.maximum(P_yx, 0.0)  # Ensure non-negative
            elif global_ref_mask is not None:
                # Use the masked similarities
                P_yx = L_masked / (tau + EPS)  # Convert back from logits if tau was applied  
                P_yx = np.maximum(P_yx, 0.0)  # Ensure non-negative
            else:
                P_yx = K_raw
            
            # Normalize P_yx to sum to (1 - P_null) for each target
            P_sum = np.sum(P_yx, axis=0)
            P_yx = P_yx * (1 - P_null) / (P_sum + EPS)
            
        else:
            # Standard logit-based approach with softmax
            # Add the Null Threshold (which competes with the masked logits)
            null_row = np.full((1, K_raw.shape[1]), tau * C_null)
            L_aug = np.vstack([L_masked, null_row])
            
            # Softmax handles -inf correctly (exp(-inf) = 0)
            P_aug = softmax(L_aug, axis=0)
            
            P_yx = P_aug[:-1, :]
            P_null = P_aug[-1, :]
        
        # Debug: Check actual effective usage after softmax
        if use_topk_masking or global_ref_mask is not None:
            # Count non-zero weights per target
            non_zero_weights = np.sum(P_yx > 1e-10, axis=0)  # Count weights > epsilon
            avg_non_zero = np.mean(non_zero_weights)
        
        # Entropy-Based Confidence
        H_norm = safe_entropy(P_yx) / np.log(P_yx.shape[0] + EPS)
        omega_k_eff = np.mean((1 - H_norm) * (1 - P_null))
        
        # Virtual Reference Projection
        w_y = (1 - P_null)
        norm_y = np.sum(w_y) + EPS
        w_y_norm = w_y / norm_y
        
        w_x = P_yx @ w_y_norm
        norm_x = np.sum(w_x) + EPS
        w_x_norm = w_x / norm_x
        
        # Gene-Specific Parameter Estimation (Vectorized - No Loop)
        # PACE v2.3: Rank Navigation, Raw Correction
        if self.hyperparams.use_v23_correction:
            # PHASE 1: NAVIGATION using RANKS (shape-based matching)
            # Generate rank space for navigation (invariant to scaling/shifts)
            X_k_ranks = np.argsort(np.argsort(X_k, axis=0), axis=0).astype(float)
            Y_k_ranks = np.argsort(np.argsort(Y_k_nav, axis=0), axis=0).astype(float)  # Use navigation data for ranks
            
            # Calculate similarity using ranks (Spearman-like correlation)
            if self.hyperparams.similarity_metric == 'centered_cosine':
                # Center the ranks
                X_k_ranks_c = X_k_ranks - np.mean(X_k_ranks, axis=0, keepdims=True)
                Y_k_ranks_c = Y_k_ranks - np.mean(Y_k_ranks, axis=0, keepdims=True)
                
                # Normalize
                X_k_ranks_norm = X_k_ranks_c / (np.linalg.norm(X_k_ranks_c, axis=0, keepdims=True) + EPS)
                Y_k_ranks_norm = Y_k_ranks_c / (np.linalg.norm(Y_k_ranks_c, axis=0, keepdims=True) + EPS)
                
                # Rank-based similarity matrix
                K_rank = X_k_ranks_norm.T @ Y_k_ranks_norm
            else:
                # Standard cosine on ranks
                X_k_ranks_norm = X_k_ranks / (np.linalg.norm(X_k_ranks, axis=0, keepdims=True) + EPS)
                Y_k_ranks_norm = Y_k_ranks / (np.linalg.norm(Y_k_ranks, axis=0, keepdims=True) + EPS)
                K_rank = X_k_ranks_norm.T @ Y_k_ranks_norm
            
            # Use rank-based similarity for P_yx calculation
            L_rank = self.hyperparams.tau * K_rank
            
            # Apply Top-K masking on rank similarities
            if self.hyperparams.use_hard_gating:
                N_ref = L_rank.shape[0]
                k_neighbors = max(5, int(N_ref * self.hyperparams.hard_gating_ratio))
                k_cutoff_idx = N_ref - k_neighbors
                
                mask = np.full_like(L_rank, -np.inf)
                for col in range(L_rank.shape[1]):
                    cutoff_val = np.partition(L_rank[:, col], k_cutoff_idx)[k_cutoff_idx]
                    keep_mask = L_rank[:, col] >= cutoff_val
                    mask[keep_mask, col] = L_rank[keep_mask, col]
                L_rank = mask
            
            # Softmax with null option
            null_row = np.full((1, L_rank.shape[1]), self.hyperparams.tau * C_null)
            L_aug = np.vstack([L_rank, null_row])
            P_aug = softmax(L_aug, axis=0)
            
            P_yx_rank = P_aug[:-1, :]
            P_null = P_aug[-1, :]
            
            # Calculate weights using rank-based P_yx
            w_y = (1 - P_null)
            w_y_norm = w_y / (np.sum(w_y) + EPS)
            w_x = P_yx_rank @ w_y_norm
            norm_x = np.sum(w_x) + EPS
            w_x_norm = w_x / norm_x
            
            # PHASE 2: CORRECTION using RAW VALUES (magnitude-preserving)
            # Calculate correction parameters using original intensity data
            mu_x_raw = np.sum(X_k * w_x_norm, axis=1)  # Raw means for correction
            var_x_raw = np.sum(w_x_norm * (X_k - mu_x_raw[:, None])**2, axis=1)
            
            mu_y_raw = np.sum(Y_k * w_y_norm, axis=1)  # Raw means for correction
            var_y_raw = np.sum(w_y_norm * (Y_k - mu_y_raw[:, None])**2, axis=1)
            
            # PHASE 3: SAFETY - Verify navigation correlation
            # Virtual reference using navigation data
            X_virtual_raw = X_k @ P_yx_rank
            
            # Correlation gate on navigation data (unscaled)
            mu_virt_raw = np.sum(X_virtual_raw * w_y_norm, axis=1)
            mu_obs_raw = np.sum(Y_k_nav * w_y_norm, axis=1)  # Use navigation data for verification
            
            X_virt_c = X_virtual_raw - mu_virt_raw[:, None]
            Y_obs_c = Y_k_nav - mu_obs_raw[:, None]  # Use navigation data for verification
            
            cov = np.sum(w_y_norm * X_virt_c * Y_obs_c, axis=1)
            v_x = np.sum(w_y_norm * X_virt_c**2, axis=1)
            v_y = np.sum(w_y_norm * Y_obs_c**2, axis=1)
            
            denom = np.sqrt(v_x * v_y) + EPS
            rho_g = cov / denom  # Raw correlation for trust
            
            trust = np.maximum(0, rho_g)**2
            alpha_raw = (np.sqrt(var_x_raw) + EPS) / (np.sqrt(var_y_raw) + EPS)
            
            # Use raw means for final correction
            mu_x = mu_x_raw
            mu_y = mu_y_raw
            
        else:
            # Original PACE v2.2 logic (correction on same data as navigation)
            X_virtual = X_k @ P_yx
            
            # Weighted Statistics (same data for navigation and correction)
            mu_x = np.sum(X_k * w_x_norm, axis=1)
            var_x = np.sum(w_x_norm * (X_k - mu_x[:, None])**2, axis=1)
            
            mu_y = np.sum(Y_k * w_y_norm, axis=1)
            var_y = np.sum(w_y_norm * (Y_k - mu_y[:, None])**2, axis=1)
            
            # Weighted correlation (same as before)
            mu_virt = np.sum(X_virtual * w_y_norm, axis=1)
            mu_obs = np.sum(Y_k * w_y_norm, axis=1)
            
            X_virt_c = X_virtual - mu_virt[:, None]
            Y_obs_c = Y_k - mu_obs[:, None]
            
            cov = np.sum(w_y_norm * X_virt_c * Y_obs_c, axis=1)
            v_x = np.sum(w_y_norm * X_virt_c**2, axis=1)
            v_y = np.sum(w_y_norm * Y_obs_c**2, axis=1)
            
            denom = np.sqrt(v_x * v_y) + EPS
            rho_g = cov / denom
            
            trust = np.maximum(0, rho_g)**2
            alpha_raw = (np.sqrt(var_x) + EPS) / (np.sqrt(var_y) + EPS)
        
        # Common logic for both versions
        # Enhanced Variance Trap - Distrust large changes in either direction
        deviation = np.abs(np.log(alpha_raw + EPS))
        large_deviation_mask = deviation > np.log(self.hyperparams.lambda_damp)
        trust[large_deviation_mask] /= self.hyperparams.lambda_damp
        
        # Final Proposals (alpha/beta calculated on original data in both cases)
        alpha_proposals = 1.0 + trust * (alpha_raw - 1.0)
        beta_proposals = mu_x - alpha_proposals * mu_y
        
        # Calculate average correlation and alpha for this pathway
        avg_correlation = np.mean(rho_g)
        avg_alpha = np.mean(alpha_proposals)
        avg_trust = np.mean(trust)
        
        # Correlation-Alpha relationship metrics
        high_corr_mask = rho_g > 0.8  # High correlation threshold
        low_corr_mask = rho_g < 0.4   # Low correlation threshold
        
        high_alpha_high_corr = np.mean(alpha_proposals[high_corr_mask]) if np.any(high_corr_mask) else 0.0
        high_alpha_low_corr = np.mean(alpha_proposals[low_corr_mask]) if np.any(low_corr_mask) else 0.0
        
        # Quality metric: correlation should predict alpha strength
        corr_alpha_correlation = np.corrcoef(rho_g, alpha_proposals)[0, 1] if len(rho_g) > 1 else 0.0
        
        # Pathway coherence analysis (eigengene correlation)
        # Calculate how well genes in this pathway move together
        try:
            from sklearn.decomposition import PCA
            
            # Center the pathway data
            X_k_centered = X_k - np.mean(X_k, axis=1, keepdims=True)
            
            if X_k_centered.shape[0] > 1 and X_k_centered.shape[1] > 1:
                # PCA on samples x genes to get pathway eigengene
                pca = PCA(n_components=1)
                pca.fit(X_k_centered.T)
                eigengene = pca.transform(X_k_centered.T).flatten()
                
                # Flip sign to align with mean expression
                if np.mean(eigengene * np.mean(X_k_centered, axis=0)) < 0:
                    eigengene = -eigengene
                
                # Calculate gene-eigengene correlations
                gene_eigengene_corrs = []
                for i in range(X_k_centered.shape[0]):
                    gene_vec = X_k_centered[i, :]
                    if np.std(gene_vec) > 1e-8 and np.std(eigengene) > 1e-8:
                        corr = np.corrcoef(gene_vec, eigengene)[0, 1]
                        gene_eigengene_corrs.append(corr)
                
                pathway_coherence = np.mean(gene_eigengene_corrs) if gene_eigengene_corrs else 0.0
                pathway_pc1_var = pca.explained_variance_ratio_[0]
                
                # Apply coherence weighting to parameter proposals
                if self.coherence_scores and pathway_name and pathway_name in self.coherence_scores:
                    # Get coherence weights for genes in this pathway
                    coherence_weights = np.ones(len(pathway_indices))
                    for i, gene_idx in enumerate(pathway_indices):
                        if hasattr(self, 'common_genes') and gene_idx < len(self.common_genes):
                            gene_name = self.common_genes[gene_idx]
                            if gene_name in self.coherence_scores[pathway_name]:
                                coherence_weights[i] = self.coherence_scores[pathway_name][gene_name]
                    
                    # Apply coherence weighting: ω_k_eff * ρ_coherence
                    # Higher coherence = more stable = higher weight
                    coherence_factor = np.mean(coherence_weights)
                    omega_k_eff *= coherence_factor
                
            else:
                pathway_coherence = 0.0
                pathway_pc1_var = 0.0
                
        except ImportError:
            # sklearn not available, skip coherence analysis
            pathway_coherence = 0.0
            pathway_pc1_var = 0.0
        except Exception:
            # Any other error in coherence calculation
            pathway_coherence = 0.0
            pathway_pc1_var = 0.0
        
        # --- ENHANCED DIAGNOSTIC CALCULATION ---
        def calc_gini(w):
            if len(w) <= 1: return 0.0
            w_sorted = np.sort(w)
            n = len(w)
            return 1.0 - 2.0 * np.sum((np.cumsum(w_sorted) - w_sorted/2.0) / np.sum(w_sorted)) / n

        def calc_neff(w):
            # Effective sample size: 1 / sum(w^2)
            return 1.0 / (np.sum(w**2) + EPS)
        
        metrics = {
            'gini_x': calc_gini(w_x_norm),
            'gini_y': calc_gini(w_y_norm),
            'neff_x': calc_neff(w_x_norm),
            'neff_y': calc_neff(w_y_norm),
            # Correlation-Correction relationship metrics
            'avg_correlation': avg_correlation,
            'avg_alpha': avg_alpha,
            'avg_trust': avg_trust,
            'high_corr_alpha': high_alpha_high_corr,
            'low_corr_alpha': high_alpha_low_corr,
            'corr_alpha_corr': corr_alpha_correlation,
            'pct_high_corr': np.mean(high_corr_mask) * 100,
            'pct_low_corr': np.mean(low_corr_mask) * 100,
            # Pathway coherence metrics
            'pathway_coherence': pathway_coherence,
            'pathway_pc1_var': pathway_pc1_var,
            # Top-K masking metrics
            'topk_masking_enabled': use_topk_masking,
            'reference_usage_pct': 100 * calc_neff(w_x_norm) / len(w_x_norm) if len(w_x_norm) > 0 else 0.0
        }
        
        return omega_k_eff, alpha_proposals, beta_proposals, K_raw.flatten(), metrics
    
    def iterative_estimation_loop_with_diagnostics(self, X_prime: np.ndarray, Y_prime: np.ndarray, 
                                                  pathway_indices_list: List[np.ndarray],
                                                  pathway_names_list: List[str],
                                                  alpha_glob: np.ndarray, beta_glob: np.ndarray,
                                                  N: int, M: int, Y_raw: np.ndarray = None) -> Tuple[np.ndarray, np.ndarray, List[Dict]]:
        """
        Step 4: The Iterative Estimation Loop with diagnostic data collection and coherence weighting.
        For PACE v3.2, includes iterative consensus mechanism.
        """
        G = X_prime.shape[0]
        C_null = 0.0
        diagnostic_data = []
        
        # PACE v3.2: Initialize consensus tracking
        if self.hyperparams.use_v32_iterative_consensus:
            consensus_S = 1.0  # Start with no global scaling
            pathway_votes_history = []
        
        for t in range(self.hyperparams.max_iter):
            # Calculate global reference mask once per iteration
            global_ref_mask = self.calculate_global_reference_mask(X_prime, Y_prime)
            
            sigma_alpha = alpha_glob * self.hyperparams.w_prior
            sigma_beta = beta_glob * self.hyperparams.w_prior
            sigma_w = np.full(G, self.hyperparams.w_prior)
            
            K_hist = []
            
            # Diagnostic Accumulators
            acc_metrics = {
                'gini_x': [], 'gini_y': [], 'neff_x': [], 'neff_y': [],
                'avg_correlation': [], 'avg_alpha': [], 'avg_trust': [],
                'high_corr_alpha': [], 'low_corr_alpha': [], 'corr_alpha_corr': [],
                'pct_high_corr': [], 'pct_low_corr': [],
                'pathway_coherence': [], 'pathway_pc1_var': [],
                'topk_masking_enabled': [], 'reference_usage_pct': []
            }
            pathway_variances = []
            informative_pathways = 0
            
            # PACE v3.2: Collect pathway votes for consensus
            if self.hyperparams.use_v32_iterative_consensus:
                pathway_alpha_votes = []
                pathway_beta_votes = []
            
            for pathway_indices, pathway_name in zip(pathway_indices_list, pathway_names_list):
                if len(pathway_indices) < 5: continue
                
                # Unpack metrics dict with pathway name for coherence weighting and Top-K masking
                omega_k_eff, alpha_props, beta_props, k_raw_vals, p_metrics = self.pathway_execution(
                    X_prime, Y_prime, pathway_indices, C_null, pathway_name, 
                    use_topk_masking=self.hyperparams.use_hard_gating,
                    global_ref_mask=global_ref_mask,
                    Y_raw=Y_raw
                )
                
                K_hist.extend(k_raw_vals)
                pathway_variances.append(np.var(k_raw_vals))
                
                if omega_k_eff > 0.01:
                    informative_pathways += 1
                    # Collect metrics only for informative pathways
                    for k, v in p_metrics.items():
                        acc_metrics[k].append(v)
                    
                    # PACE v3.2: Collect pathway votes
                    if self.hyperparams.use_v32_iterative_consensus:
                        # Calculate pathway-specific global shift estimates
                        pathway_alpha_mean = np.mean(alpha_props)
                        pathway_beta_mean = np.mean(beta_props)
                        pathway_alpha_votes.append(pathway_alpha_mean)
                        pathway_beta_votes.append(pathway_beta_mean)
                
                for i, gene_idx in enumerate(pathway_indices):
                    sigma_alpha[gene_idx] += alpha_props[i] * omega_k_eff
                    sigma_beta[gene_idx] += beta_props[i] * omega_k_eff
                    sigma_w[gene_idx] += omega_k_eff
            
            # PACE v3.2: Consensus Decision Making
            if self.hyperparams.use_v32_iterative_consensus and len(pathway_alpha_votes) > 3:
                # Calculate consensus statistics
                alpha_votes = np.array(pathway_alpha_votes)
                beta_votes = np.array(pathway_beta_votes)
                
                alpha_mean = np.mean(alpha_votes)
                alpha_std = np.std(alpha_votes)
                beta_mean = np.mean(beta_votes)
                beta_std = np.std(beta_votes)
                
                # Consensus threshold: Low variance = Technical artifact, High variance = Biological
                alpha_consensus = alpha_std < self.hyperparams.consensus_threshold
                beta_consensus = beta_std < self.hyperparams.consensus_threshold
                
                pathway_votes_history.append({
                    'iteration': t + 1,
                    'alpha_votes': alpha_votes.tolist(),
                    'beta_votes': beta_votes.tolist(),
                    'alpha_mean': alpha_mean,
                    'alpha_std': alpha_std,
                    'beta_mean': beta_mean,
                    'beta_std': beta_std,
                    'alpha_consensus': alpha_consensus,
                    'beta_consensus': beta_consensus,
                    'n_votes': len(alpha_votes)
                })
                
                # Apply consensus if pathways agree
                if alpha_consensus and abs(alpha_mean - 1.0) > 0.1:  # Significant shift detected
                    # Update global scaling for next iteration
                    new_S = consensus_S * alpha_mean
                    print(f"  🗳️  Consensus detected: {len(alpha_votes)} pathways vote α={alpha_mean:.3f}±{alpha_std:.3f}")
                    print(f"  🔧 Updating global S: {consensus_S:.3f} → {new_S:.3f}")
                    
                    # Re-scale Y_prime for next iteration
                    if Y_raw is not None:
                        Y_prime = np.arcsinh(np.sinh(Y_raw) * new_S)
                    consensus_S = new_S
                else:
                    print(f"  🤔 No consensus: {len(alpha_votes)} pathways, α_std={alpha_std:.3f} (threshold={self.hyperparams.consensus_threshold})")
            
            if K_hist:
                C_new = np.percentile(K_hist, 5)
                C_null = (1 - self.hyperparams.eta) * C_null + self.hyperparams.eta * C_new
            
            # --- CALCULATE AGGREGATE DIAGNOSTICS ---
            def safe_mean(lst): return np.mean(lst) if lst else 0.0
            
            avg_gini_x = safe_mean(acc_metrics['gini_x']) # REFERENCE Sparsity (The Goal)
            avg_neff_x = safe_mean(acc_metrics['neff_x'])
            avg_reference_usage_pct = safe_mean(acc_metrics['reference_usage_pct'])
            
            # Override reference usage calculation if using global masking
            if global_ref_mask is not None:
                global_avg_refs_per_target = np.mean(np.sum(global_ref_mask, axis=0))
                avg_reference_usage_pct = 100 * global_avg_refs_per_target / N
            
            # Calculate ratios
            ratio_x = avg_neff_x / N if N > 0 else 0.0
            
            diagnostic_record = {
                'iter': t + 1,
                'tau': self.hyperparams.tau,
                'mean_pathway_gini': avg_gini_x, 
                'neff_x': avg_neff_x,
                'neff_x_ratio': ratio_x,
                'gini_x': avg_gini_x,
                'N': N, 'M': M,
                'n_pathways_total': len(pathway_indices_list),
                'n_pathways_informative': informative_pathways,
                'c_null': C_null,
                'reference_usage_pct': avg_reference_usage_pct
            }
            
            # Add consensus diagnostics for v3.2
            if self.hyperparams.use_v32_iterative_consensus:
                diagnostic_record['consensus_S'] = consensus_S
                if pathway_votes_history:
                    latest_vote = pathway_votes_history[-1]
                    diagnostic_record['alpha_consensus'] = latest_vote['alpha_consensus']
                    diagnostic_record['alpha_std'] = latest_vote['alpha_std']
                    diagnostic_record['n_pathway_votes'] = latest_vote['n_votes']
            
            diagnostic_data.append(diagnostic_record)
            
            # Minimal diagnostic printing - only final iteration
            if t == self.hyperparams.max_iter - 1:
                if self.hyperparams.use_v32_iterative_consensus:
                    print(f"  Final: Ref_Usage={ratio_x:.1%}, Gini_X={avg_gini_x:.3f}, Pathways={informative_pathways}/{len(pathway_indices_list)}, Consensus_S={consensus_S:.3f}")
                else:
                    print(f"  Final: Ref_Usage={ratio_x:.1%}, Gini_X={avg_gini_x:.3f}, Pathways={informative_pathways}/{len(pathway_indices_list)}")

        # Final aggregation
        EPS = 1e-8
        total_weight = np.sum(sigma_w)
        if total_weight <= len(sigma_w) * self.hyperparams.w_prior + EPS:
            alpha_final = np.ones(G) + np.random.normal(0, 0.01, G)
            beta_final = np.random.normal(0, 0.01, G)
        else:
            alpha_final = sigma_alpha / (sigma_w + EPS)
            beta_final = sigma_beta / (sigma_w + EPS)
        
        # Store consensus history for v3.2
        if self.hyperparams.use_v32_iterative_consensus:
            self.consensus_history = pathway_votes_history
        
        return alpha_final, beta_final, diagnostic_data
    
    def calculate_gene_pathway_coherence(self, X_data: np.ndarray, gene_indices: np.ndarray) -> Dict[str, Dict[str, float]]:
        """
        Calculate gene-pathway coherence (null stability prior) for all pathways.
        This measures how well each gene correlates with its pathway eigengene.
        
        Args:
            X_data: Reference data matrix (genes x samples)
            gene_indices: Gene names corresponding to rows
            
        Returns:
            Dict mapping pathway_name -> {gene_name: coherence_score}
        """
        
        # Map genes to matrix indices
        gene_map = {g: i for i, g in enumerate(gene_indices)}
        coherence_scores = {}
        
        for pathway_name, gene_list in self.pathway_dict.items():
            # Get pathway indices
            pathway_indices = [gene_map[g] for g in gene_list if g in gene_map]
            
            if len(pathway_indices) < 5:  # Skip small pathways
                continue
                
            # Extract pathway data
            X_pathway = X_data[pathway_indices, :]
            
            try:
                from sklearn.decomposition import PCA
                
                # Center the data
                X_centered = X_pathway - np.mean(X_pathway, axis=1, keepdims=True)
                
                if X_centered.shape[0] > 1 and X_centered.shape[1] > 1:
                    # Calculate pathway eigengene (PC1)
                    pca = PCA(n_components=1)
                    pca.fit(X_centered.T)
                    eigengene = pca.transform(X_centered.T).flatten()
                    
                    # Flip sign to align with mean expression
                    if np.mean(eigengene * np.mean(X_centered, axis=0)) < 0:
                        eigengene = -eigengene
                    
                    # Calculate gene-eigengene correlations
                    pathway_coherence = {}
                    for i, gene_idx in enumerate(pathway_indices):
                        gene_name = gene_indices[gene_idx]
                        gene_vec = X_centered[i, :]
                        
                        if np.std(gene_vec) > 1e-8 and np.std(eigengene) > 1e-8:
                            corr = np.corrcoef(gene_vec, eigengene)[0, 1]
                            # Convert correlation to stability score (0-1 range)
                            # High correlation = high stability = high weight
                            stability = max(0.0, corr)  # Only positive correlations count
                            pathway_coherence[gene_name] = stability
                        else:
                            pathway_coherence[gene_name] = 0.0
                    
                    coherence_scores[pathway_name] = pathway_coherence
                    
            except ImportError:
                # sklearn not available, use uniform weights
                pathway_coherence = {gene_indices[i]: 1.0 for i in pathway_indices}
                coherence_scores[pathway_name] = pathway_coherence
            except Exception:
                # Any error, use uniform weights
                pathway_coherence = {gene_indices[i]: 1.0 for i in pathway_indices}
                coherence_scores[pathway_name] = pathway_coherence
        
        return coherence_scores
    
    def save_diagnostic_data(self, diagnostic_data: List[Dict], output_path: str) -> None:
        """Save diagnostic data to CSV file with robust error handling."""
        import os
        
        # Ensure output directory exists - FAIL FAST if cannot create
        output_dir = os.path.dirname(output_path)
        try:
            os.makedirs(output_dir, exist_ok=True)
        except Exception as e:
            raise RuntimeError(f"CRITICAL: Cannot create diagnostic output directory {output_dir}: {str(e)}")
        
        # Verify directory is writable
        if not os.access(output_dir, os.W_OK):
            raise RuntimeError(f"CRITICAL: Diagnostic output directory is not writable: {output_dir}")
        
        if pd is not None:
            # Use pandas if available
            try:
                df = pd.DataFrame(diagnostic_data)
                df.to_csv(output_path, index=False)
            except Exception as e:
                raise RuntimeError(f"CRITICAL: Failed to save diagnostic data with pandas: {str(e)}")
        else:
            # Fallback to manual CSV writing
            try:
                if not diagnostic_data:
                    with open(output_path, 'w') as f:
                        f.write("# No diagnostic data available\n")
                    print(f"WARNING: No diagnostic data to save, created empty file: {output_path}")
                    return
                
                # Get all unique keys from all records
                all_keys = set()
                for record in diagnostic_data:
                    all_keys.update(record.keys())
                all_keys = sorted(all_keys)
                
                # Write CSV manually
                with open(output_path, 'w') as f:
                    # Header
                    f.write(','.join(all_keys) + '\n')
                    
                    # Data rows
                    for record in diagnostic_data:
                        values = [str(record.get(key, '')) for key in all_keys]
                        f.write(','.join(values) + '\n')
            except Exception as e:
                raise RuntimeError(f"CRITICAL: Failed to save diagnostic data manually: {str(e)}")
        
        # Verify file was created and has content
        if not os.path.exists(output_path):
            raise RuntimeError(f"CRITICAL: Diagnostic file was not created: {output_path}")
        
        file_size = os.path.getsize(output_path)
        if file_size == 0:
            raise RuntimeError(f"CRITICAL: Diagnostic file is empty: {output_path}")
        
        print(f"✓ Diagnostic data successfully saved: {output_path} ({file_size} bytes, {len(diagnostic_data)} records)")
    
    def extrapolate_unique_genes(self, alpha_common: np.ndarray, beta_common: np.ndarray,
                               mu_Y_common: np.ndarray, mu_Y_unique: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """
        Extrapolate parameters for unique genes using lookup table interpolation.
        """
        if len(mu_Y_unique) == 0:
            return np.array([]), np.array([])
        
        if len(mu_Y_common) == 0:
            return np.ones_like(mu_Y_unique), np.zeros_like(mu_Y_unique)
        
        # Sort common genes by intensity for interpolation
        sort_idx = np.argsort(mu_Y_common)
        mu_sorted = mu_Y_common[sort_idx]
        alpha_sorted = alpha_common[sort_idx]
        beta_sorted = beta_common[sort_idx]
        
        # Interpolate parameters for unique genes
        alpha_unique = np.interp(mu_Y_unique, mu_sorted, alpha_sorted)
        beta_unique = np.interp(mu_Y_unique, mu_sorted, beta_sorted)
        
        return alpha_unique, beta_unique
    
    def align(self, ref_data: BatchData, target_data: BatchData, 
             diagnostic_output_path: Optional[str] = None) -> Tuple[np.ndarray, Dict]:
        """
        Main PACE v2.2 alignment function with corrected gene mapping and diagnostic output.
        
        Args:
            ref_data: Reference batch data
            target_data: Target batch data  
            diagnostic_output_path: Optional path to save diagnostic CSV file
        """
        print("Starting PACE v2.2 alignment...")
        
        # Initialize diagnostics collection
        diagnostic_data = []
        
        # Step 1: Find common and unique genes using intersect1d with return_indices
        common_genes, x_indices, y_indices = np.intersect1d(
            ref_data.gene_indices, target_data.gene_indices, return_indices=True
        )
        
        # Extract data for common genes
        X_common = ref_data.data[x_indices]
        Y_common = target_data.data[y_indices]
        
        # Find unique genes in target
        unique_mask = ~np.isin(target_data.gene_indices, common_genes)
        Y_unique = target_data.data[unique_mask]
        unique_genes = target_data.gene_indices[unique_mask]
        
        # Step 2: Adaptive preprocessing
        if self.hyperparams.use_v35_affine_anchors:
            # PACE v3.5: Global affine anchors (returns affine parameters)
            X_prime, Y_prime, Y_raw, affine_params = self.adaptive_preprocessing(X_common, Y_common)
            S = affine_params  # Store affine params for metadata
            alpha_glob, beta_glob = affine_params
            print(f"Global affine transform: α={alpha_glob:.4f}, β={beta_glob:.4f}")
        elif self.hyperparams.use_v34_housekeeping_anchors:
            # PACE v3.4: Housekeeping anchors (returns Y_raw for navigation safety)
            X_prime, Y_prime, Y_raw, S = self.adaptive_preprocessing(X_common, Y_common)
        elif self.hyperparams.use_v32_iterative_consensus:
            # PACE v3.2: Iterative consensus (starts with S=1.0)
            X_prime, Y_prime, S = self.adaptive_preprocessing(X_common, Y_common)
            Y_raw = np.arcsinh(Y_common)  # Keep unscaled version for consensus updates
        elif self.hyperparams.use_v31_pure_local:
            # PACE v3.1: Pure local estimation (no global scaling)
            X_prime, Y_prime, S = self.adaptive_preprocessing(X_common, Y_common)
            Y_raw = Y_prime  # Same as Y_prime since no scaling applied
        elif self.hyperparams.use_v30_activity:
            X_prime, Y_prime, Y_raw, S = self.adaptive_preprocessing(X_common, Y_common)
        else:
            preprocessing_result = self.adaptive_preprocessing(X_common, Y_common)
            if len(preprocessing_result) == 4:  # v2.4 case
                X_prime, Y_prime, Y_raw, S = preprocessing_result
            else:  # Standard case
                X_prime, Y_prime, S = preprocessing_result
                # PACE v2.4: Create unscaled target data for navigation
                if self.hyperparams.use_v24_navigation:
                    Y_raw = np.arcsinh(Y_common)  # Unscaled for navigation
                else:
                    Y_raw = Y_prime  # Use scaled data for both navigation and correction (v2.2/v2.3 behavior)
        
        if not self.hyperparams.use_v35_affine_anchors:
            print(f"Global scale factor S: {S:.4f}")
        
        # Step 3: Initialize with ComBat baseline
        alpha_glob, beta_glob = self.combat_baseline.compute_baseline(X_prime, Y_prime)
        
        # Step 4: Map pathways to gene indices and calculate coherence scores
        common_gene_map = {g: i for i, g in enumerate(common_genes)}
        pathway_indices_list = []
        pathway_names_list = []
        
        for pathway_name, gene_list in self.pathway_dict.items():
            indices = [common_gene_map[g] for g in gene_list if g in common_gene_map]
            if len(indices) >= 5:  # Minimum pathway size 5 (as in corrected version)
                pathway_indices_list.append(np.array(indices))
                pathway_names_list.append(pathway_name)
        
        # Calculate coherence scores for null stability prior
        if len(pathway_indices_list) > 0:
            self.common_genes = common_genes  # Store for coherence weighting
            self.coherence_scores = self.calculate_gene_pathway_coherence(X_prime, common_genes)
        
        if len(pathway_indices_list) == 0:
            alpha_final = alpha_glob
            beta_final = beta_glob
            
            # Record diagnostic data for baseline case
            diagnostic_data.append({
                'iter': 0,
                'tau': self.hyperparams.tau,
                'neff_x': ref_data.data.shape[1],  # All samples equally weighted
                'neff_y': target_data.data.shape[1],
                'neff_x_ratio': 1.0,  # Uniform weights
                'neff_y_ratio': 1.0,
                'gini_x': 0.0,  # No sparsity
                'gini_y': 0.0,
                'N': ref_data.data.shape[1],
                'M': target_data.data.shape[1],
                'n_pathways_total': 0,
                'n_pathways_informative': 0,
                'pathway_variance_mean': 0.0,
                'pathway_variance_std': 0.0,
                'pathway_variance_max': 0.0
            })
        else:
            # Step 5: Iterative estimation loop with diagnostic collection and coherence weighting
            alpha_final, beta_final, diagnostic_data = self.iterative_estimation_loop_with_diagnostics(
                X_prime, Y_prime, pathway_indices_list, pathway_names_list, alpha_glob, beta_glob,
                ref_data.data.shape[1], target_data.data.shape[1], Y_raw
            )
        
        # Step 6: Handle unique genes
        if len(unique_genes) > 0:
            # Calculate mean intensity on prime data for interpolation
            mu_Y_common = np.mean(Y_prime, axis=1)
            # Transform unique to prime scale for interpolation x-axis
            Y_unique_prime = np.arcsinh(Y_unique * S)
            mu_Y_unique = np.mean(Y_unique_prime, axis=1)
            
            alpha_unique, beta_unique = self.extrapolate_unique_genes(
                alpha_final, beta_final, mu_Y_common, mu_Y_unique
            )
        else:
            alpha_unique, beta_unique = np.array([]), np.array([])
        
        # Step 7: Assemble final parameters with optimized gene mapping
        full_alpha = np.ones(len(target_data.gene_indices))
        full_beta = np.zeros(len(target_data.gene_indices))
        
        # Create a map for ALL target indices once: O(N)
        target_gene_map = {g: i for i, g in enumerate(target_data.gene_indices)}
        
        # 1. Fill Common Genes (Vectorized map lookup)
        for gene, val_alpha, val_beta in zip(common_genes, alpha_final, beta_final):
            if gene in target_gene_map:
                idx = target_gene_map[gene]
                full_alpha[idx] = val_alpha
                full_beta[idx] = val_beta
        
        # 2. Fill Unique Genes
        if len(unique_genes) > 0:
            for gene, val_alpha, val_beta in zip(unique_genes, alpha_unique, beta_unique):
                if gene in target_gene_map:
                    idx = target_gene_map[gene]
                    full_alpha[idx] = val_alpha
                    full_beta[idx] = val_beta
        
        # Step 8: Apply clamping
        full_alpha = np.clip(full_alpha, 0.1, 10.0)
        
        # Step 9: Final application
        if self.hyperparams.use_v35_affine_anchors:
            # For v3.5, S contains affine parameters (alpha, beta)
            alpha_glob, beta_glob = S
            # Apply the same affine transform that was used in preprocessing
            # Y_corrected = (Y - beta) * alpha, then apply local corrections
            Y_affine_corrected = (target_data.data - beta_glob) * alpha_glob
            Y_prime_corrected = full_alpha[:, np.newaxis] * np.arcsinh(Y_affine_corrected) + full_beta[:, np.newaxis]
            Y_final = np.sinh(Y_prime_corrected)
        else:
            # Standard application for other versions
            Y_prime_corrected = full_alpha[:, np.newaxis] * np.arcsinh(target_data.data * S) + full_beta[:, np.newaxis]
            Y_final = np.sinh(Y_prime_corrected)
        
        # Step 10: Save diagnostic data if path provided
        if diagnostic_output_path:
            self.save_diagnostic_data(diagnostic_data, diagnostic_output_path)
        
        # Return results
        metadata = {
            "S": S,
            "n_common_genes": len(common_genes),
            "n_unique_genes": len(unique_genes),
            "n_pathways": len(pathway_indices_list),
            "clamps": (0.1, 10.0),
            "hyperparams": self.hyperparams,
            "diagnostics": diagnostic_data
        }
        
        print("PACE v2.2 alignment completed.")
        return Y_final, metadata

# ==========================================
# Legacy PACE class for backward compatibility
# ==========================================

class PACE(PACE_v22):
    """Legacy PACE class that inherits from PACE_v22 for backward compatibility."""
    pass