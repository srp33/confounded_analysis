Biological Context for Bulk Transcriptomic Batch Correction

1. The Nature of the Data

Microarray vs. RNA-seq

While the method is platform-agnostic, the input data ($X$ and $Y$) typically originates from two distinct technologies with unique noise profiles:

Microarray (Intensity-based): Measures continuous fluorescence intensity.

Characteristics: Limited dynamic range (saturation at high expression, background noise at low expression). Data is continuous but often heteroscedastic (variance increases with mean).

Artifacts: Probe-specific binding affinities, hybridization efficiency.

Bulk RNA-seq (Count-based): Measures discrete read counts mapping to gene features.

Characteristics: Theoretically infinite dynamic range, but subject to "dropout" in low-expressed genes. Data is strictly non-negative and often modeled as Negative Binomial.

Artifacts: Library depth (sequencing depth), GC-content bias, gene length bias.

The "Gene Manifold" Assumption

Genes do not behave independently. They operate in tightly regulated co-expression networks (pathways).

Manifold Hypothesis: High-dimensional gene expression data lies on a lower-dimensional manifold defined by biological constraints.

Implication for Correction: A valid batch correction method must respect these correlations. If Gene A and Gene B are strongly correlated in the Reference, they should remain correlated in the Corrected Target unless biological differences dictate otherwise.

2. The Nature of Batch Effects

The "Orthogonality" Fallacy

A common, yet dangerous, assumption in standard batch correction (e.g., ComBat, PCA regression) is that Batch Effects are orthogonal to Biological Signal.

The Assumption: The vector of batch variation is mathematically perpendicular to the vector of biological variation.

The Reality: In meta-analyses and observational studies, batch and biology are often confounded.

Example: Batch 1 contains mostly "Control" samples (Healthy). Batch 2 contains mostly "Case" samples (Disease).

Consequence: "Blind" methods that regress out the batch effect will inevitably regress out the biological signal, stripping the "Disease" signature because it looks like a batch difference.

Topological Distortion vs. Linear Shifts

Batch effects manifest in two primary ways:

Linear Shift/Scale (Global): The entire distribution is shifted up/down (mean offset) or stretched/compressed (variance difference). This is common due to normalization differences (e.g., TPM vs FPKM).

Non-Linear Distortion (Topology): The relationships between samples are warped.

The Danger: Standard methods (like Quantile Normalization) force distributions to be identical. If Batch 1 is bimodal (two distinct cell types) and Batch 2 is unimodal (one cell type), forcing them to match destroys the biological truth of Batch 1.

3. Expectations for the Agent

The AI must understand that:

Preservation is Priority: We prefer under-correction (leaving some noise) to over-correction (erasing biology).

Heterogeneity is Real: Samples are rarely pure. They are mixtures of cell types. A "sample" is a vector summation of its constituent cell profiles.

Anchors must be Semantic: We cannot blindly match Sample A to Sample B based on raw correlation alone, because raw correlation includes the batch effect. We must find "Semantic Anchors" (shared biological processes) to define similarity.