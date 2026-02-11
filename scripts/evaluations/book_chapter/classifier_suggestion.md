Yes — and the fact that **plain logistic is performing very poorly is actually a strong diagnostic signal**, not just a nuisance.

Given *your* regime (small (n), large (p), correlated features, batch-corrected expression data), there **is one classifier that clearly deserves to replace logistic** rather than just sit alongside it.

---

## Short answer (the one to add)

> **Add Regularized Quadratic Discriminant Analysis (RDA)**
> *(and feel free to drop plain logistic entirely)*

If I had to add **exactly one** classifier to *replace* logistic in your list, it would be **RDA**.

---

## Why logistic is failing (this matters)

Logistic regression fails badly when:

* Class separation is **not linearly separable in the original feature space**
* Signal is **distributed across correlated genes**
* Marginal effects are weak but **joint structure matters**
* Batch correction introduces subtle nonlinear distortions

Elastic net fixes *some* of this, but it is still:

* Linear
* Discriminative
* Coefficient-based

So logistic performing poorly is expected and *informative*.

---

## Why RDA is the right replacement

### What RDA does

RDA interpolates smoothly between:

| Model       | What it assumes           |
| ----------- | ------------------------- |
| Naive Bayes | Independent features      |
| LDA         | Shared covariance         |
| QDA         | Class-specific covariance |

By shrinking class covariances toward a common (or diagonal) structure, RDA:

* Explicitly models **correlations**
* Allows **class-dependent geometry**
* Remains stable in small-n settings

This is exactly the failure mode of logistic.

---

## Why RDA is **not redundant** with shrinkage LDA

| Shrinkage LDA     | RDA                                     |
| ----------------- | --------------------------------------- |
| Linear boundary   | **Quadratic boundary**                  |
| Shared covariance | **Partially class-specific covariance** |
| Strong bias       | Lower bias                              |
| Very stable       | Slightly more variance                  |

If:

* Shrinkage LDA is good but not great → RDA often improves it
* Elastic net beats LDA sometimes → RDA often closes the gap

---

## Practical recommendation in your context

### Keep:

* Elastic Net
* SVM
* RF
* XGBoost (when stable)
* Shrinkage LDA

### Replace:

* **Logistic → RDA**

Logistic adds *no new inductive bias* beyond elastic net in your regime.

---

## How to add RDA cleanly in R

Use **`klaR::rda`**:

```r
library(klaR)

rda_fit <- klaR::rda(
  x = X_train,
  grouping = y_train,
  gamma = 0.5,  # covariance shrinkage (0 = LDA, 1 = QDA)
  lambda = 0.5  # diagonal shrinkage
)

pred <- predict(rda_fit, X_test)

probs <- pred$posterior[, 2]
```

Start with:

```r
gamma = 0.3–0.7
lambda = 0.3–0.7
```

You do *not* need heavy tuning.

---


## Final verdict

> **Logistic performing poorly is telling you that linear discriminative models are insufficient.**
> **RDA directly addresses that failure mode and is the cleanest replacement.**

5. Concrete RDA hyperparameters for your case

Given 
𝑛
=
300
,
𝑝
=
15,000
n=300,p=15,000:

gamma  = 0.3  # how class-specific covariance is
lambda = 0.6  # how diagonal the covariance is


Interpretation:

Mostly shared covariance

Heavy shrinkage