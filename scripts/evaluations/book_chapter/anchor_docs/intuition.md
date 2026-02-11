Project Intuition: The Philosophy of Semantic Batch Correction

1. The Core Problem: The "Frankenstein" Dataset

Traditional batch correction methods (like MNN or CCA) operate on mathematical proximity in high-dimensional space.

The Flaw: In high dimensions, distance is unreliable. A "T-cell" in Batch A might be mathematically closer to a "NK-cell" in Batch B due to technical noise than it is to a "T-cell" in Batch B.

The Result: "Blind" alignment stitches these incorrect samples together, creating a "Frankenstein" dataset where cell identities are blurred or swapped.

2. The Solution: Semantic Anchoring

This project rejects "Blind" mathematical alignment in favor of Semantic Alignment.

The Intuition: Instead of asking "Are these two samples numerically close?", we ask "Do these two samples are activate the same biological programs?"

The Mechanism: We use Gene Sets (Pathways) as the fundamental unit of comparison.

If Sample $i$ (Reference) and Sample $j$ (Target) both show high activity in Glycolysis, Apoptosis, and T-cell Receptor Signaling, they are biologically similar, regardless of the global batch shift.

3. The Power of the Ensemble (Consensus)

A single pathway (e.g., "Cell Cycle") might be corrupted by batch effects (e.g., if the batch effect mimics proliferation signals).

The "Wisdom of Crowds": By aggregating votes from hundreds of distinct pathways, random noise cancels out, while true biological signal (which is consistent across multiple functional axes) is amplified.

Robustness: This makes the method resistant to outliers and "confounded" batch effects.

4. Why Linear Transformation (Moment Matching)?

We explicitly choose a Linear Transformation ($y' = \alpha y + \beta$) over non-linear warping (like Quantile Normalization or warping functions).

Preserving Bimodality (The "Valley" Problem)

Scenario: Consider a gene (e.g., ESR1 in breast cancer) that is bimodal: some samples are negative, some are positive. There is a deep "valley" between the two peaks.

The Risk: Aggressive non-linear methods (like Quantile Normalization) often "fill in" this valley to force the distribution to match a Gaussian reference.

Our Approach: Linear scaling ($\alpha$) stretches the distribution, and shifting ($\beta$) moves it, but the topology (the valley) is preserved. If the data implies two distinct populations, our method keeps them distinct.

5. Goals

Correct Technical Bias: Remove the global shift and scale differences between platforms.

Preserve Biological Heterogeneity: Maintain the distinct clusters and non-normal distributions inherent in biology.

Explainability: The correction weights are derived from known biology (pathways), making the correction traceable, unlike opaque "black box" neural network autoencoders.