Diagnostic & Validation Protocol: The CPW-BC "Stress Test" Suite

1. Philosophy of Testing: "Mixing is Easy, Preservation is Hard"

Standard batch correction metrics (like kBET or LISI) focus on "mixing" — ensuring that Batch 1 and Batch 2 overlap.

The Trap: A method that sets all values to zero achieves perfect mixing (kBET = 0) but destroys all biology.

The Mandate: This protocol prioritizes Biological Preservation over Batch Removal. We assume under-correction is safer than over-correction.

2. Test 1: The "Valley" Preservation Test (Topology)

Objective: Verify that the linear correction ($\alpha y + \beta$) preserves the bimodal structure of key marker genes, rather than forcing them into a Gaussian shape (as ComBat or Quantile Normalization might).

Procedure:

Identify Sentinel Genes: Select top 50 genes with high bimodal coefficients (Ashman's D or Dip Test > 0.05) in the Reference batch ($X$).

Measure Topology:

Calculate the Kurtosis and Skewness of these genes in $X$.

Calculate the same metrics in the Corrected Target ($Y'$).

Metric: Topology_Retention_Score = Correlation(Kurtosis_X, Kurtosis_Y')

Pass: Score > 0.8. (The "shape" of the distribution is preserved).

Fail: Score < 0.5. (The correction has "normalized" the biology away).

3. Test 2: The "Frankenstein" Monitor (Semantic Integrity)

Objective: Ensure that the "Consensus" is finding legitimate biological anchors, not just forcing matches between unrelated samples to minimize distance.

Diagnostic Plot: Trust vs. Deviation

X-axis (Trust): The Semantic Similarity Score ($W_{ij}$) derived from the pathway ensemble.

Y-axis (Deviation): The magnitude of correction applied ($|\alpha_{local} - 1|$).

Interpretation:

Healthy Triangle: You should see high deviation only when Trust is high. (i.e., "We are correcting this sample heavily, but only because we are 99% sure it's a match").

The "Hallucination" Zone: High deviation at Low/Moderate Trust.

Flag: If samples with $W_{ij} < 0.5$ are receiving $\alpha$ corrections deviating > 20% from unity, the algorithm is "guessing" and creating Frankenstein profiles.

4. Test 3: The Semi-Synthetic "Ground Truth" Test

Since real data lacks ground truth, we must simulate it to verify the math.

Procedure:

Construct Ground Truth ($Y_{true}$): Take the Reference Batch ($X$). Select a random subset of samples to serve as the "True Target".

Condition A: The Recovery Test (Injected Artifacts)

Inject Batch Effect:

Linear: Add a global shift (+2.0) and scale ($\times 1.5$).

Non-Linear: Apply a "Saturation" curve to high-expression genes ($y = \tanh(x)$).

Dropout: Randomly set 20% of low values to zero.

Run CPW-BC: Attempt to recover $Y_{true}$ from the corrupted version.

Metrics:

RMSE: Euclidean distance between $Y'$ and $Y_{true}$.

Parameter Recovery: Can the model recover the injected Linear Scale (1.5) and Shift (2.0)?

Condition B: The Null Control (No Artifacts)

Setup: Use $Y_{true}$ directly as the input Target without modification.

Run CPW-BC: Correct $Y_{true}$ against $X$.

Expectation: The method should detect that distributions match and apply identity mapping ($\alpha=1, \beta=0$).

Metrics:

Silence Score: Mean absolute deviation of $\alpha$ from 1.0 and $\beta$ from 0.0.

Pass: Deviation < 0.01. (If the method corrects data that needs no correction, it is over-fitting/hallucinating).

5. Test 4: The Orphan Gene Consistency Check

Objective: Verify that "Orphan Genes" (genes not in any pathway) are not being left behind or corrected incoherently.

Procedure:

Partition Genes: Split genes into Set A (Pathway Genes) and Set B (Orphans).

Compare Distributions: Plot the distribution of the scaling factors ($\alpha$) for Set A vs. Set B.

Expectation:

If the distributions are distinct (e.g., Mean $\alpha_A = 1.2$, Mean $\alpha_B = 1.0$), the method is failing to propagate the correction to the rest of the transcriptome.

Fix: This diagnostic triggers the need for "Gene Manifold Neighbor Borrowing" (using correlations in Reference to assign orphan genes to pathway proxies).

6. Test 5: Robustness to Outliers (The "Explosion" Check)

Objective: Test the fragility of the Moment Matching equation ($\alpha = \sigma_X / \sigma_Y$).

Procedure:

Injection: Take a single gene in the Target batch. Change one sample's value to $100 \times$ the max.

Run Correction: Measure $\alpha$ for that gene.

Fail Condition: If $\alpha$ drops by > 50% due to a single sample, the method is too brittle.

Remediation: If this test fails, replace Standard Deviation ($\sigma$) with Median Absolute Deviation (MAD) in the math:

$$\alpha = \frac{\text{MAD}(X)}{\text{MAD}(Y)}$$

7. Test 6: The Compositional Imbalance Challenge (Biology Preservation)

Objective: Verify that the algorithm distinguishes between "Batch Effect" (global technical shift) and "Compositional Difference" (biological shift). Standard quantile normalization fails here because it forces the global distribution of Batch B to look like Batch A, effectively crushing the dominant cell type in B if it's rare in A.

Procedure:

Source Data: Use a dataset with at least two distinct biological clusters (e.g., Cell Type A and Cell Type B).

Create Skewed Batches:

Reference ($X$): Subsample to create a mix of 90% Type A / 10% Type B.

Target ($Y$): Subsample to create a mix of 10% Type A / 90% Type B.

Inject Batch Effect: Add a large global shift (+3.0) to Target ($Y$).

Run Correction: Correct $Y$ to match $X$.

Evaluation:

Failure Mode (Simpson's Paradox): If the corrected $Y$ shifts the mean of Type B (dominant in $Y$) to match the global mean of $X$ (dominated by Type A), the method has failed. It treated biology as batch.

Success: The corrected $Y$ Type B samples should align with the rare $X$ Type B samples, not the abundant $X$ Type A samples.

Metric:

Cluster Alignment Error: Calculate the Euclidean distance between the centroid of $X_{TypeB}$ and $Y'_{TypeB}$.

Pass: Distance < 0.2 (standard deviations).

8. Summary Scorecard

The agent should output a final report card:

Metric

Description

Target

Valley Score

Correlation of Kurtosis (Topology)

> 0.8

Frankenstein Index

% of Low-Trust samples with High-Correction

< 5%

Linear Recovery

Accuracy of retrieving injected gain/offset

> 95%

Silence Score

Deviation from identity when no batch effect exists

< 0.01

Orphan Alignment

Divergence between Pathway & Orphan correction

< 0.1

Outlier Stability

Change in $\alpha$ per unit of outlier sigma

< 0.1

Compositional Safety

Error in aligning rare-to-dominant clusters

< 0.2