  
The reference CSV has 267 samples (all 3 training studies including USA itself), but we're computing batch labels for only 197 (excluding USA). The `adjust_target_data_SA_UK.R` script uses all `n` studies as reference without excluding the test study. Fix is to derive batch labels by matching sample names from the CSV against `dat_lst`:

Editedcheck_batch_class_split.R+14-5

Also fix the class balance summary which uses `tidyr` without loading it:

Editedcheck_batch_class_split.R+1-1

```
pixi run Rscript scripts/check_batch_class_split.R
```

Shouldn't we be reducing the effects of confounding by using the class labels for adjustment? Curious.

Yes, and that's exactly the right intuition — and also exactly why the result is so counterintuitive. The key is that supervision helps the _training_ correction but creates a harder problem at _projection time_.

Here's the mechanism step by step:

**What supervised ComBat is doing correctly:** in the training data, if batch 1 (GSE37250_SA) is mostly Active TB and batch 2 (India) is mostly Control, unsupervised ComBat would see the gene expression difference between batches and "correct" it — but it's actually removing biological signal, not batch effects. Supervised ComBat uses the class labels to disentangle the two: "the difference between batches _after accounting for class composition_ is the batch effect." This is genuinely correct, and the training space it creates is better — class-separated, batch-mixed.

**Where it goes wrong:** the training space is now _calibrated to the training batches' specific class compositions_. When the test data is projected into this space (step 2, unsupervised ComBat, no test labels), ComBat aligns the test batch's marginal distribution to the training's marginal distribution. But the training marginal was constructed assuming specific class proportions across batches. If the test batch has different class composition, the alignment shifts test samples systematically into the wrong region. The more "class-purified" the training space (higher PC1 variance), the worse this misalignment is.

**Why KNN is uniquely hurt:** a linear classifier can still recover — it learns weights on the whole feature space and can compensate for a systematic shift by adjusting its intercept. KNN has no such mechanism. It just looks up nearest neighbors in raw feature space, and if the test batch has been shifted into the wrong-class cluster, every test sample's neighbors are wrong.

**The deeper irony:** supervised ComBat was designed for the use case of correcting a single dataset for downstream analysis — not for a train/test setting where test data must be projected without labels. Using it in a train/test pipeline essentially punishes you for having confounded training data: the stronger the confounding, the more aggressively supervised ComBat corrects it, the more class-calibrated the training space becomes, and the more vulnerable the unsupervised test projection is.

So yes — supervised ComBat is doing exactly what it's supposed to do for the training data. The failure is in the projection step, which has no equivalent "supervised" mechanism. The fix would either be to use test labels during projection (leakage) or to use a method that doesn't require projecting test data into a labeled training space at all.

Hmm, there's a clear difference between /home/phr23/confounded_analysis/scripts/evaluations/book_chapter/outputs/diagnostics/combat_sup_knn/class_within_batch_combat_sup.png and /home/phr23/confounded_analysis/scripts/evaluations/book_chapter/outputs/diagnostics/combat_sup_knn/class_within_batch_combat.png

Read 2 files

This is the smoking gun. Look at what combat_sup has done to the India and USA panels: each batch has been collapsed into a **single tight cluster**, separated by class along PC1 — India Control on the left (~~-20), India Active TB on the right (~~+20), USA Control on the left (~~-20), USA Active TB on the right (~~+20).

Compare to combat: all three batches are spread out and intermixed along PC1, with no clear class separation per batch.

What supervised ComBat has done is dramatic: it has taken each batch and **squeezed all within-class variance to near zero**, then placed the two class clusters at opposite ends of PC1. PC1 now explains 50% of variance (up from 25.3% under combat) because the entire space has been reorganised around the class label.

This is the mechanism behind the KNN failure in full:

1. Supervised ComBat finds that within each training batch, the class label explains most variance, so it "preserves" that by amplifying it and removing everything else — collapsing each batch into two tight class-specific points in space.
2. The test data (USA) is then projected via unsupervised ComBat into this space. The projection has to match the test batch's marginal distribution to the training marginal. But the training marginal is now dominated by these four tight clusters (two batches × two classes), not a smooth distribution.
3. The test samples all get mapped onto those four cluster positions — and which cluster they land on is determined purely by which training batch they're most similar to in the raw space, **not by their actual class label**.

This also explains why the effect is so extreme (purity ~0.19, worse than random): the within-class variance has been so thoroughly collapsed that there's almost no information left about the true class of the test samples — they get assigned to clusters based on batch similarity alone.

If you were to modify Combat, (see /home/phr23/BioPreserve/basis_workflow/scripts/ComBat.R), what changes would you make to prevent this problem? Or is it unpreventable? In my mind, we should be able to better adjust the training data by using the labels.

Let me read the ComBat implementation to ground any suggestion in the actual code.

ReadComBat.R

Having read it, my honest assessment: the training correction in this code is _fine_ — the failure is almost entirely in how test data gets projected afterward, and parts of it are fixable while one part is fundamental.

Let me separate what's happening mechanically, because it points to the specific changes.

## Where the damage actually happens

The supervised training fit is correct. Look at lines 172–178: `stand.mean` is built as `grand.mean + class_effect` (the covariate columns of `B.hat`), so the class signal is moved into `stand.mean` and _protected_ from the batch estimation. Then `var.pooled` (lines 156–170) is the **residual within-class, within-batch variance** — it's small precisely because the class signal has been pulled out into the mean. That's all working as intended.

The catastrophe is in your step 2 — the separate, unsupervised `ComBat(mod=NULL)` call that projects the test batch. Two distinct pathologies:

1. **Mean collapse.** With `mod=NULL`, `stand.mean` for the test batch is just the class-_marginalized_ grand mean. Test Active-TB and test Control both get pulled toward that single mean, which sits _between_ the two training class clusters. So test samples land in no-man's-land or biased toward whichever class dominates training — the wrong-neighbor problem.
2. **Variance over-shrinkage.** Line 288 divides the test batch by `sqrt(delta.star)`, and `delta.star` is EB-shrunk toward the training pooled variance — which is the _tiny within-class residual_. The test data still carries its class signal as part of its total variance, so matching it to that tiny residual collapses exactly the signal KNN needs. This is what produced the worse-than-random 0.19 purity.

## Changes, ranked

**(a) `mean.only=TRUE` for the test batch only — cheap, low-risk.** Skip the `delta.star` scaling at line 288 for the projected batch. This directly stops pathology #2. The test data keeps its own spread instead of being crushed to the training residual scale. This alone should pull KNN back from worse-than-random to roughly neutral.

**(b) Project against class-agnostic reference statistics.** Instead of a fresh unsupervised ComBat, reuse the training fit's `grand.mean` and `var.pooled` and only remove the test batch's _mean_ offset from the grand mean. You're aligning the test batch to the class-marginalized centroid and pooled scale — never trying to reproduce the bimodal class structure you can't see. This addresses #1 without inventing labels.

**(c) Semi-supervised projection (the principled one).** Pseudo-label the test set — a quick classifier trained on the corrected training data, or even k-means with k=2 — and pass those as the `mod` covariate when correcting the test batch. This restores the symmetry that's currently broken: training is covariate-protected, test should be too. Soft pseudo-labels (class probabilities) avoid hard leakage and degrade gracefully when wrong.

## Is it fundamentally preventable?

Partly. The **information asymmetry is real and unfixable**: you cannot place test samples into a space whose geometry is _defined_ by labels without using labels. So supervised-ComBat-then-blind-projection will never _beat_ plain unsupervised ComBat for a distance method like KNN — your intuition that labels should help is right for the training data, but the benefit can't survive a label-free projection.

What _is_ fixable is the **active harm**. Worse-than-random purity isn't the asymmetry talking — it's variance/mean shrinkage actively anti-aligning the test data. Fix (a) and (b) and you should get back to "neutral, no worse than unsupervised." Fix (c) and you can recover some of the benefit, bounded by how good your pseudo-labels are.

If you want, I can prototype (a)+(b) as a `combat_sup_v2` adjuster in `adjust_target_data_SA_UK.R` and re-run the KNN-purity diagnostic to confirm it climbs back above 0.5.

Let's project against class-agnostic reference statistics.

Let me read the current combat_sup implementation to see exactly how to modify the projection step.

Readadjust_target_data_SA_UK.R

Now I understand the current structure. Note line 467 actually passes `combined_group` (test labels) into the projection — but regardless, it aligns the test batch to the frozen bimodal training reference, which is what scatters the test samples.

I'll implement option (b) as a **new** adjuster `combat_sup_ca` (class-agnostic) so we can A/B compare against the original `combat_sup` rather than overwrite it. Let me add it to the valid adjusters and implement the method.

Editedadjust_target_data_SA_UK.R+29-1

Now generate the adjusted data for the n3/testUSA scenario so we can compare:

RanGenerate class-agnostic combat_sup adjusted data for n3 testUSA

Data generated. Now let me write a focused comparison that runs KNN purity for all three methods (combat, combat_sup, combat_sup_ca) on this scenario:

Created a file, ran a command

combat_sup_ca made it _worse_ (0.174 vs 0.334), not better. That's an important negative result. Let me look at the PCA to understand why:

Readcombat_sup_ca_pca.png

This is a clear and informative negative result: **combat_sup_ca did not help — it's actually slightly worse** (0.174 vs 0.334 at k=5). The PCA reveals why, and it means my diagnosis of _where_ the problem lives was wrong.

Look at the `combat_sup_ca` panel. The training data still has its bimodal structure — Control cluster on the left (PC1 ≈ −20), Active TB on the right (PC1 ≈ +20). But the test triangles are placed wrong: many **red triangles (test Active TB) sit on the left, near the training Control cluster**, and blues are scattered through the middle. The test's class signal is anti-aligned with training's class axis, which is exactly why purity falls below random.

Here's the corrected understanding. The problem is **not** in the projection statistics (mean/variance matching) — fixing those can't help. The problem is **upstream, in the training reference geometry itself**:

Your training is batch-confounded — GSE37250_SA is mostly Active TB, India is mostly Control. When supervised ComBat estimates the class effect `B.hat` from this data, that estimate is **confounded with the SA-vs-India batch difference**. The "class axis" it builds (and amplifies) is really a _class-tangled-with-batch_ axis. Within each training batch the classes do separate along it (which is why the earlier `class_within_batch` plot looked clean), but the _absolute position_ of the class clusters in the corrected space encodes batch identity.

A new batch (USA) doesn't share SA's or India's batch signature, so when its samples are scored on this contaminated axis, they land on the wrong side — regardless of how carefully we compute the projection. My class-agnostic affine map faithfully preserved the test's own internal structure, then dropped it into a reference frame whose axes don't mean what we assumed.

This reframes the whole finding:

- **Plain `combat` works (0.88 purity)** precisely _because_ it doesn't build a contaminated class axis — it only removes batch location/scale and leaves the biology where it naturally falls.
- **`combat_sup` fails** because it amplifies a batch-confounded class axis, and the amplification makes test misplacement catastrophic.
- **No projection-side fix can recover this**, because the damage is baked into the reference geometry before projection happens.

The honest conclusion: for cross-study classification with **batch-confounded training data**, supervised batch correction is contraindicated for distance-based classifiers — the supervision imports the training confound into the feature geometry. The actionable rule would be to gate supervised ComBat on training-batch class balance: it's only safe when each training batch is reasonably class-balanced, so `B.hat` isn't confounded.

Want me to verify this directly — e.g., measure the correlation between the training class axis (`B.hat` covariate direction) and the SA-vs-India batch-mean difference vector? If they're highly aligned, that's quantitative proof the class axis is batch-contaminated.

Yes please

Let me write a script that quantifies the contamination directly — measuring how aligned the supervised class axis is with the batch-contrast axis, and crucially whether it generalizes to the test study's actual class direction.

Created a file, ran a command

This **refutes my contamination hypothesis** — and corrects a false premise I had. Two things jump out:

1. **The training pair isn't confounded at all.** SA is 46 Active / 48 Control, India is 53 Active / 50 Control — both essentially balanced. I was wrong that "SA is mostly Active, India mostly Control." So the confounding story was built on a faulty assumption.
2. **The supervised class axis is clean and generalizes well.** `cos(class_adj, batch_contr) = −0.17` (barely any batch contamination), and critically `cos(class_adj, class_USA) = +0.79` — the supervised class direction aligns strongly with USA's _actual_ class direction in raw gene space. The naive axis generalizes about equally well (+0.76).

So the class signal points the **right way** and transfers to the test set. That means the KNN failure is **not** a direction problem. It must be a **positioning** problem: even with correctly-oriented class clusters, if the test cloud is _shifted or scaled_ along the axes relative to training, every test sample can still land nearest to wrong-class training points.

This is a sharper, more answerable question. The direct test: project the corrected train and test data onto the class axis (computed from corrected training) and look at where the four groups — train-Active, train-Control, test-Active, test-Control — actually fall.

Created a file, ran a command

This is the answer. Look at the **class ordering** table — it's the smoking gun:

|method|Train (Active − Control)|Test (Active − Control)|
|---|---|---|
|combat|**+24.5**|**+24.2**|
|combat_sup|**+20.1**|**−5.65**|
|combat_sup_ca|+20.1|−10.8|

Under `combat`, Active TB scores higher than Control on the class axis in _both_ train and test — consistent, so KNN works. Under `combat_sup`, the test ordering **flips sign**: test Active TB scores _lower_ than test Control, the opposite of training.

The absolute positions make the KNN failure mechanical. For combat_sup along the class axis:

- Train Control ≈ −15.7, Train Active ≈ +4.5
- Test Active ≈ −8.4 → nearest training cluster is **Control** (−15.7) ✗
- Test Control ≈ −2.8 → nearest training cluster is **Active** (+4.5) ✗

Both test classes map to the wrong training cluster — purity below random, exactly as observed.

The crucial part: we already proved the _raw_ test class direction is correct (`cos(class_adj, class_USA) = +0.79`). So **the inversion is introduced by combat_sup's adjustment itself**, specifically the step-2 projection — not by any pre-existing property of the data.

**Why step 2 inverts it.** ComBat's correction is per-gene affine: `corrected_g = a_g + b_g · raw_g`. The projection onto the class axis `w` becomes `Σ w_g·b_g·raw_g` — the per-gene **variance rescaling `b_g` reweights the genes that define the class axis**. When the amplified class structure dominates training variance, step-2's `delta` rescaling differentially down/up-weights those genes enough to flip the net projection. This is the same variance-shrinkage pathology, now shown to actually invert the class signal.

That also explains why `combat_sup_ca` didn't help — I targeted the reference _statistics_, but the inversion comes from per-gene _rescaling_ during projection, which my affine map inherited.

**The fix this points to:** make step-2 location-only (`mean.only=TRUE`). With `b_g = 1` for all genes, the projection becomes `const + w·raw_test` — the class ordering is preserved up to a shift, so it can't invert. This is concrete and directly testable.

Want me to implement `combat_sup` with a mean-only step-2 projection and re-run the class-axis + KNN-purity check to confirm the inversion disappears?

Give it a try—but I expect that method will fail, if for example, the two training sets are microarray, and the test set is rna-seq.

You're absolutely right to flag that. Mean-only step-2 preserves ordering but abandons scale correction entirely — and for cross-platform data, scale differences are real and large. What we actually need is scale correction that can't flip a direction: a **single global scale factor** rather than per-gene `delta_g`. A positive scalar applied uniformly can never invert the class axis, and it still handles the microarray-vs-RNA-seq dynamic range problem.

Read a file, edited a file

Now add the implementation after `combat_sup_ca`:

Editedadjust_target_data_SA_UK.R+40-0

Generate the adjusted data and run the comparison:

Ran 2 commands, edited 2 files

The class-axis ordering is still inverted for `combat_sup_mg` (−15.4), and KNN purity climbs to 0.489 but doesn't recover. The global scale correction helped meaningfully compared to the other variants — but the ordering is still wrong.

The inversion is confirmed to be coming from step 2's **mean-only location shift itself**, not delta. Look at what mean-only ComBat is doing to the test data: after the supervised-ComBat training, the train class centroids are at −15.7 (Control) and +4.5 (Active). The mean of that distribution is around −5 (more samples near Control). When mean-only step-2 shifts the test data to match that mean, it pulls the entire test cloud leftward — and since the test's raw class signal is weaker than training's amplified signal, both test classes end up to the left of the training Active centroid, with Active landing nearer Control.

In short: **the training reference frame is so class-polarized** by supervised ComBat (class mean separation = 20 units) that any projection of the unlabeled test batch — regardless of whether you use mean-only or per-gene scale — will place the test cloud in a region that misaligns with at least one class. The geometry is simply incompatible with a label-free projection step.

This is the fundamental limit. The only way to project correctly into a class-polarized training space is to know which class each test sample belongs to. Any label-free projection strategy is working against the geometry that supervised ComBat deliberately created.

So to directly answer your original question: supervised ComBat **is** using the labels correctly to improve training, but that very improvement makes the training space impossible to project into without labels. The problem is architectural — you cannot separate "good supervised correction of training" from "bad unlabeled projection of test" as long as the projection target is the supervised-corrected training distribution. The fix would require either accepting that KNN doesn't pair with supervised ComBat, or using a projection target that wasn't class-polarized — i.e., projecting test against the _unsupervised_-corrected training distribution for step 2 and then applying only the mean-shift difference from supervised correction. That's genuinely a new variant worth testing.