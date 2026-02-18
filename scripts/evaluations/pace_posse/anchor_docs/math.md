Mathematical Specification: Consensus Pathway-Weighted Moment Matching

1. Global Definitions

$X \in \mathbb{R}^{G \times N}$: Reference Matrix (Genes $\times$ Samples).

$Y \in \mathbb{R}^{G \times M}$: Target Matrix (to be corrected).

$\mathcal{S} = \{S_1, S_2, \dots, S_K\}$: A collection of $K$ biological pathways, where each $S_k$ is a set of gene indices.

2. Preprocessing & Variance Stabilization

Data is first transformed to stabilize variance, typically using the inverse hyperbolic sine, which linearizes log-normal data while handling zeros gracefully.

$$X_{stab} = \text{asinh}(X), \quad Y_{stab} = \text{asinh}(Y)$$

3. Pathway-Level Subspace Projection

For each pathway $S_k$, we project the global data into a local subspace. This filters out noise from irrelevant genes.

$$X_k = X[S_k, :], \quad Y_k = Y[S_k, :]$$

4. Local Similarity Kernels (The "Votes")

Within each subspace $k$, we compute a similarity matrix $C_k \in \mathbb{R}^{N \times M}$ representing the affinity between every Reference sample and every Target sample based only on that pathway.

$$C_k(i, j) = \text{CosineSimilarity}(X_k^{(i)}, Y_k^{(j)})$$

Note: High-Contrast Filtering (Softmax with Temperature $\tau$) is often applied here to sharpen the distinction between true matches and noise.

5. Consensus Weight Aggregation

We aggregate the local kernels into a global consensus weight matrix $W$. This relies on the "Ensemble Hypothesis": individual pathways may be noisy, but the consensus across hundreds of pathways is robust.

$$W_{total} = \sum_{k=1}^{K} C_k$$

To ensure robust matching, we apply Mutual Nearest Neighbor (MNN) logic (or "Bi-stochastic" normalization) to $W_{total}$. This prevents "Hubness" (where one Reference sample matches everything).

$$P_{ij} = \sqrt{ P(i|j) \cdot P(j|i) }$$

Where $P(i|j)$ and $P(j|i)$ are column-wise and row-wise normalized versions of $W_{total}$.

6. Parameter Estimation: Weighted Moment Matching

Instead of using unstable regression (OLS/WLS), we use Closed-Form Weighted Moments. This effectively aligns the first two moments (Mean and Variance) of the Target distribution to the Reference distribution, weighted by biological similarity.

For each gene $g$:

A. Compute Weighted Statistics

We compute the weighted mean ($\mu$) and weighted standard deviation ($\sigma$) of the Reference ($X$) and Target ($Y$) using the Consensus Weights $W$.

$$\mu_{X,g} = \frac{\sum_i W_{i} \cdot X_{g,i}}{\sum_i W_{i}}, \quad \sigma_{X,g} = \sqrt{\frac{\sum_i W_{i} (X_{g,i} - \mu_{X,g})^2}{\sum_i W_{i}}}$$

(Similarly for $Y$, using weights derived from the column sums of $W$).

B. Solve for Scaling ($\alpha$) and Shift ($\beta$)

We seek a linear transformation $y' = \alpha y + \beta$ such that the moments match.

Scale (Gain) Parameter:


$$\alpha_g = \frac{\sigma_{X,g}}{\sigma_{Y,g}}$$

Shift (Offset) Parameter:


$$\beta_g = \mu_{X,g} - (\alpha_g \cdot \mu_{Y,g})$$

Why Moment Matching? It is analytically closed-form, computationally efficient, and avoids the "shrinkage" problems of ridge regression. It preserves the shape (kurtosis/skewness) of the original distribution.

7. Final Correction

The correction is applied element-wise to the Target matrix:

$$Y'_{g,j} = \alpha_g \cdot Y_{g,j} + \beta_g$$

Finally, the inverse transformation is applied to return to linear space:


$$Y_{final} = \sinh(Y')$$