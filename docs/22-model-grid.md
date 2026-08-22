# Model grid, fixed in advance

Referenced by `docs/21-preregistration.md` §6. Written before any model was
fitted. Selection happens on the **validation** period (2010–2016) only.

## Preprocessing (not tuned)

- Features standardised to zero mean and unit variance using **training-period
  statistics only**. Applying validation or test statistics would leak.
- Constant columns (zero training variance) are dropped and the count reported.
  They carry no information and would divide by zero.
- No feature selection, no dimensionality reduction, no correlation pruning. The
  point of the exercise is to give the model everything; pruning by any criterion
  is a physical judgement smuggled in through the back door.

## A. Conditional logistic regression

Stratified on the matched set: within each stratum of 1 case and 10 controls the
model predicts which row is the case, so the likelihood is a softmax over the
stratum and the intercept cancels.

| hyperparameter | grid |
|---|---|
| L2 strength λ | 1e-4, 1e-3, 1e-2, 1e-1, 1, 10, 100, 1000 |

**8 configurations.** Optimiser: L-BFGS to convergence (gradient norm < 1e-6) or
500 iterations, whichever first. The objective is convex, so the optimum does not
depend on initialisation and there is no restart to tune.

## B. Gradient-boosted trees

Ranking within strata, same softmax objective, so the two model classes are
scored on exactly the same likelihood and their information gains are comparable.

| hyperparameter | grid |
|---|---|
| trees | 100, 300, 1000 |
| max depth | 2, 3, 4 |
| learning rate | 0.01, 0.05 |

**18 configurations.** Subsample 0.8, minimum 20 rows per leaf — fixed, not tuned.

## Total

26 configurations. Whichever has the highest validation information gain becomes
the primary model, and that decision is made and recorded before the test period
is opened. The full validation table for all 26 is reported regardless of outcome,
so that the gap between best and median is visible — if the best config beats the
median by roughly the amount you would expect from noise across 26 draws, that is
itself the finding.

## What is deliberately not here

- Early stopping on validation *and* reporting validation as the result. Early
  stopping is allowed; the stopped model is then scored on test, once.
- Ensembling the top-k configurations. That is a 27th model chosen after seeing
  26 scores.
- Any per-magnitude, per-region or per-depth model. Those are subgroups, and
  §10 of the pre-registration forbids promoting them.
