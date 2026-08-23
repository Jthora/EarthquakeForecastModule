# Pre-registration: chart features as earthquake-timing predictors

**Written 2026-08-22, before any model was fitted to these features.**
Committed before the first training run so that its timestamp is checkable.

This document exists because of `docs/07-research-log.md`. Six times in this
programme an analysis produced a confident wrong answer, the worst reporting
p = 10⁻⁸⁹ where the true value was 0.70. Every one of them was found by a check
run *after* the result looked good. A design with ~9,800 features and ~12,000
events will fit anything asked of it; the only defence is to fix the question,
the metric, and the threshold in advance, and then to look exactly once.

---

## 1. The question

Given the positions and derived angular relationships of the Sun, Moon and
planets — the full astrological chart, computed but not interpreted — can a model
identify **when** an earthquake occurred, among candidate times at the same place,
better than chance, on data it has never seen?

This is deliberately not a question about mechanism. No physical account is
required, offered, or tested. If prediction is possible, prediction is the result.

**What this cannot answer.** The design conditions on an earthquake having
occurred in that cell. It therefore says nothing about *where* earthquakes happen
or about absolute rates — only about timing. A positive result would be a
component of a forecast, `λ(x,t) = μ(x) · f(features)`, not a forecast.

## 2. Data

| | |
|---|---|
| Catalogue | USGS ComCat, `data/comcat/global_m40.csv`, 488,215 events |
| Span | 1976-01-01 to 2024-12-31 (49.0 years) |
| Magnitude | M ≥ 5.5 primary; M ≥ 5.0 and M ≥ 4.0 as secondary scales |
| Declustering | Gardner–Knopoff, largest-first. 23,258 → **12,160** independent at M5.5+ |
| Cells | Equal-area, ~100 km, 50,930 total, 3,786 occupied at M5.5+ |

Declustering is applied before anything else and is not a tunable choice. The
windows are generous and discard some independent events; that costs power and
cannot manufacture signal, which is the correct direction to err.

## 3. Design

Matched case-control. Each case is an earthquake at its exact time; its controls
are **the same cell at other times**. Sharing the cell cancels tectonic setting,
station density and regional magnitude bias exactly, without modelling them.

| | primary | secondary |
|---|---|---|
| scheme | `DayOffset`, ±1..5 whole days | `Window`, uniform ±365 days |
| controls per case | 10 | 10 |
| also conditions out | local solar time, season, network trend | season, network trend |
| therefore blind to | diurnal effects; anything slower than ~5 days | anything slower than ~1 year |
| seed | 20260822 | 20260823 |

Two schemes are run because neither is trustworthy alone. Whole-day offsets hold
local solar time fixed, so the daily cycle in detection threshold — traffic,
quarry blasts, cultural noise — is identical for case and control and cannot be
mistaken for signal. That same property makes the primary design blind to any
genuine diurnal effect, and nearly blind to slow planetary configurations. The
±365-day window sees slow features but reacquires seasonal and network-trend
exposure. **A result that appears under one and not the other is an artefact of
that scheme's specific blindness, and will be reported as such.**

## 4. Features

9,816 per row, from `planetary-harmonics-core`, generated with no selection,
filtering, or physical justification:

| family | count | |
|---|---|---|
| aspects | 7,920 | cos/sin of n(λᵢ−λⱼ), 55 pairs × 24 harmonics × 3 frames |
| declination | 495 | parallel/contraparallel/difference |
| fixed points | 240 | aspects to the galactic centre and anticentre |
| resonance | 216 | CosmicCypher base-N scores |
| motion | 201 | speeds, retrograde flags, station proximity |
| lunar | 86 | synodic/anomalistic/draconic phases, distance, declination, node |
| station timing | 44 | days since and until each planet's station |
| chart shape | 27 | circular concentration, gaps, span, clusters |
| eclipse | 7 | true angular separation from syzygy, radius ratio, umbral margin |
| site-local | 580 | sidereal time, ascendant, midheaven, hour angles, altitudes |

Frames: geocentric, heliocentric, barycentric. Every feature is a cos/sin pair or
a rotation-invariant statistic, so **no result can depend on where the zodiac's
zero point is placed**. That is a property of the construction, not a check.

## 5. Split

Fixed now, by date, and never re-drawn.

| | period | independent events (M5.5+) |
|---|---|---|
| train | 1976-01-01 – 2009-12-31 | ~8,400 |
| validate | 2010-01-01 – 2016-12-31 | ~1,900 |
| **test** | 2017-01-01 – 2024-12-31 | ~1,900 |

The test period is **sealed**. It will be scored once, after the model and all
hyperparameters are frozen on validation. If it is scored a second time, that
fact will be recorded in the results and the result downgraded to exploratory.

A time split, not a random one, because a random split would leak: two events
hours apart in the same region would land on opposite sides and the model would
score well by recognising the region and the epoch, which is not forecasting.

## 6. Model

Primary: **L2-regularised conditional logistic regression**, stratified on the
matched set. This is the correct likelihood for the design — it conditions on
exactly one case per stratum and so cannot be fooled by the case rate.

Secondary: **gradient-boosted trees** ranking within strata, to allow the
interactions a linear model cannot see. The user's question is about *emergent*
capability, so a model class that can only find additive effects is not enough
on its own.

Whichever of the two has the higher **validation** information gain becomes the
primary reported model. That choice is made before the test set is opened.

Hyperparameters (L2 strength; trees, depth, learning rate) are selected on the
validation period only, over a grid fixed in `docs/22-model-grid.md` before
fitting. The number of configurations tried is recorded and reported.

## 7. Metric

**Information gain in bits per event**, on the sealed test set.

A stratum has 1 case and 10 controls, so a model that knows nothing assigns each
row 1/11 and the null log-likelihood per stratum is log₂(1/11) = −3.459 bits.
The reported quantity is

```
IG  =  mean over test strata of [ log₂ p̂(case) − log₂(1/11) ]
```

Bits per event, because it is the honest currency of a forecast: it is what a
user of the forecast gains, it is comparable across magnitude thresholds and
control counts, and unlike AUC or accuracy it cannot be inflated by the class
balance. Accuracy is not reported as a headline for exactly that reason — a
model that always guesses "control" is 91% accurate and worthless.

## 8. Success criterion, fixed in advance

Both conditions must hold.

1. **Detection.** Test-set IG > 0 with p < 0.01 under the block-shift permutation
   null (§9). One-sided; a negative IG is failure, not a finding.
2. **Meaningfulness.** Test-set IG ≥ **0.01 bits per event**. Below this the
   effect may be real and is still useless: 0.01 bits is roughly a 3% improvement
   in picking the true time out of eleven. Anything smaller will be reported as
   "consistent with zero for any practical purpose" regardless of its p-value.

**Replication.** The primary and secondary designs must agree in sign. If they do
not, the result is reported as scheme-dependent — that is, as evidence of
temporal bias rather than of signal.

**Declared failure.** If test IG ≤ 0, or p ≥ 0.01, or IG < 0.01 bits, the
conclusion recorded is that these features do not predict earthquake timing at
this magnitude and scale. That will be written up as plainly as a success would
be, in `docs/19-results.md`, and no post-hoc subgroup will be substituted for it.

## 9. Null calibration — run before the model, not after

The pipeline must demonstrate it produces nothing from nothing before it is
allowed to produce something.

- **Synthetic null.** A catalogue with times uniform in the span and positions
  uniform on the sphere, run through the entire pipeline. Required: IG
  indistinguishable from zero, and the p-value distribution over 200 repetitions
  KS-uniform. *(The sampler-level version of this is already a passing test:
  `sampling::tests::a_signal_free_catalogue_produces_no_case_control_separation`.)*
- **Label shuffle.** Real features, case/control labels permuted within stratum.
  Required: IG ≈ 0.
- **Block shift.** Real catalogue, event times shifted by a common offset per
  spatial block, preserving clustering while destroying any astronomical
  alignment. 200 draws give the null distribution against which the p-value in
  §8 is computed. Required: KS-uniform on the synthetic null.

**If block-shift is not KS-uniform, no result is reported at all.** A miscalibrated
null is how p = 10⁻⁸⁹ happened before.

## 10. Multiple comparisons

Every number below multiplies the chances of a false positive and is declared now:

- 2 designs (day-offset, window)
- 3 magnitude thresholds (5.5, 5.0, 4.0)
- 2 model classes (conditional logistic, boosted trees)

= **12 primary analyses.** The headline is the pre-specified primary alone:
M5.5+, day-offset, whichever model class won on validation. The other eleven are
reported as secondary with Holm correction applied across all twelve. No
subgroup, feature family, or time window not listed here will be promoted to a
headline result.

## 11. What would make this wrong anyway

Recorded now so it cannot be rationalised later.

- **Residual clustering.** Gardner–Knopoff is imperfect. A surviving aftershock
  sequence inside a ±5-day window is the single most likely source of a false
  positive. Mitigation: report results with and without a stricter declustering,
  and check that IG does not scale with the fraction of cases in dense clusters.
- **Catalogue time errors.** Origin times before ~1990 in remote regions can be
  seconds to minutes off. This blurs fast features and cannot create signal.
- **UT1 ≈ UTC.** Bounded at ±0.9 s, 0.004° of Earth rotation. Negligible.
- **Magnitude-threshold drift.** The network improved over 49 years. Day-offset
  matching removes this; the window design does not, at ±365 days.
- **Depth.** Not used in matching. Deep and shallow events respond differently to
  anything tidal; this is a known limitation, not a planned analysis.

---

*If the result is negative, this document is what makes the negative worth
having. If it is positive, this document is the only reason anyone should
believe it.*

---

## Amendment 1 — 2026-08-22, before any real fit

**Controls per case reduced from 10 to 4 for the fitted model.**

Reason: hardware, not results. The machine has 9 GB of RAM and ~2.5 GB of real
headroom. Holding the training split at 10 controls needs 3.6 GB, which drove the
machine into swap — 14 GB of it — and the fit made 2m43s of progress in 12m26s of
wall clock. Subsampling to 4 controls brings the training split to 1.65 GB.

Cost: a matched set with 1 case and *k* controls carries `k/(k+1)` of the
information available at *k* = ∞, so 10 → 4 gives up about 12% of the efficiency.
The null baseline changes from log₂(1/11) = −3.459 to log₂(1/5) = −2.322 bits;
information gain is defined relative to each stratum's own size, so the metric and
the 0.01-bit threshold in §8 are unaffected and remain comparable.

The controls kept are the first four per stratum in generator order, which is an
unbiased subsample and is deterministic given the dataset seed. **The 10-control
matrices are unchanged on disk** and can be refitted without regeneration if more
memory becomes available.

This amendment is written before any model has been fitted to real features. The
only quantities observed at this point are row counts, feature counts, timings,
and the null-calibration results on synthetic data.
