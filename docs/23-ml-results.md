# Chart features as earthquake-timing predictors — results

Pre-registered in `docs/21-preregistration.md` (2026-08-22, before any fit).
Model grid in `docs/22-model-grid.md`. **The sealed 2017–2024 test period has not
been opened.** Everything below is training and validation.

---

## 1. What was built

9,816 features per (epoch, site), generated with no selection or physical
filtering, from `planetary-harmonics-core`:

aspects (7,920: 55 pairs × 24 harmonics × cos/sin × 3 frames) · declination
parallels (495) · galactic-centre aspects (240) · CosmicCypher resonance (216) ·
motion and stations (245) · lunar synodic/anomalistic/draconic phases, distance,
node (86) · chart shape statistics (27) · eclipse geometry (7) · site-local
sidereal time, ascendant, midheaven, hour angles and altitudes (580).

Every feature is a cos/sin pair or a rotation-invariant statistic, so no result
can depend on where the zodiac's zero point sits. Frames are geocentric,
heliocentric and barycentric.

Validated against external standards, not just internally: GMST reproduces the
defined 18h41m50.548s at J2000; the Greenwich-noon Sun-to-MC offset comes out at
0.771°, which is the equation of time for 1 January; the lunar node regresses at
−19.34°/yr; 893 planetary stations are found across 49 years where orbital
periods predict ~891.

## 2. Catalogue

| | |
|---|---|
| ComCat 1976–2024 | 488,215 events |
| declustered (Gardner–Knopoff) | **154,302** independent at M4.0+, 12,160 at M5.5+ |
| removed as dependent | 68.4% at M4.0+, 47.7% at M5.5+ |
| equal-area cells | 50,930 of ~100 km; 11,090 occupied at M4.0+ |

## 3. Three bugs that silently suppressed signal

Recorded because each one produced a *plausible* null, and two of them would
have been reported as results.

**Conditioning.** Standardisation was folded into the coefficients — algebraically
correct, since the mean term cancels inside the stratum softmax. It left L-BFGS
facing raw columns whose variances run from 0.5 (a cosine) to 4×10⁸ (lunar
distance in km). At λ=1e-4, where 9,781 features against 8,414 strata should
overfit to nearly log₂5 bits, it reported convergence at 0.0005. The synthetic
null calibration missed it entirely because those features were all N(0,1) — *a
calibration set better-conditioned than the real data can only catch logical
faults, not numerical ones.*

**Float32 centring.** An outer-planet aspect changes by less than float32 epsilon
across a 5-day stratum, so its centred residual is pure quantisation noise and
z = U/√V explodes. Permutation median max|z| was **21.5**. In float64 it is 3.8.

**Non-exchangeable referents.** See §5 — the largest of the three.

## 4. Model A — conditional logistic, M5.5+, day-offset design

| λ | train | validate |
|---|---|---|
| 1e-4 | +0.253 | −4.422 |
| 1e-2 | +0.254 | −4.327 |
| 1 | +0.243 | −1.844 |
| 100 | +0.184 | −0.281 |
| 1000 | +0.146 | **−0.108** |

Every λ negative on held-out data, improving monotonically toward λ→∞ (all
coefficients zero, IG exactly 0). Best achievable: predict nothing.

**This result is uninformative, and the check that shows why is the important
part.** Planting a synthetic effect of known size into the *real* feature matrix
and refitting:

| planted effect | validation IG | |
|---|---|---|
| 0.02 (2% per SD) | −0.217 | invisible |
| 0.05 | −0.178 | invisible |
| 0.10 | −0.205 | invisible |
| 0.20 | −0.137 | invisible |
| **0.50 (50% per SD)** | −0.043 | **invisible** |

A design that cannot see a 50% modulation cannot be used to rule out a 3% one.
The null from Model A is a statement about the model, not about astrology. With
9,781 features and 8,414 strata, L2 must shrink so hard that real signal goes
with the noise.

## 5. The design flaw

The score scan's dispersion diagnostic, on M4.0+:

```
permutation z sd: 1.0014 pooled; per draw 1.001 ± 0.031
observed z sd    0.628, at permutation quantile 0.005
```

Twelve standard deviations narrow. Cause: with controls at *t* ± k days the
controls straddle the case, so for any smoothly-varying feature **the case sits
at the temporal centroid of its own stratum** and is closer to the stratum mean
than a random member. The rows are not exchangeable; every test assuming they
are is mis-specified.

The direction is conservative, so the nulls below stand — but power was being
discarded for nothing, and a null from a mis-calibrated test is not reportable.

Replaced by the **time-stratified** referent design (`Amendment 2`): referents
are every day in the case's own calendar month differing from it by a whole
multiple of 7 days. The set is fixed by the calendar rather than by the case, so
exchangeability holds by construction. Seven days also pins day of week, so the
weekly cycle in cultural noise cannot masquerade as signal.

## 6. Score scan — the powerful test

Under the matched null the case is equally likely to be any row of its stratum,
so `z_k = U_k/√V_k` is standard normal exactly, with nothing fitted. Full design
power behind each feature individually. Multiplicity handled by permuting which
row is the case, which preserves the correlation between cos and sin of
neighbouring harmonics rather than assuming independence.

| dataset | strata | max \|z\| | null median | 95th pct | family-wise p |
|---|---|---|---|---|---|
| M5.5+ day-offset | 10,154 | 3.514 | 3.791 | 4.398 | 0.90 |
| M4.0+ day-offset | 119,620 | 3.474 | 3.885 | 4.571 | 0.92 |
| M4.0+ time-stratified | — | *running* | | | |

In both completed scans the largest observed statistic is **below the median of
its own null**. No feature separates earthquakes from their referents.

## 7. Model B — gradient-boosted trees

Built (conditional stratified-softmax objective, so information gains are
directly comparable to Model A). On M5.5+ day-offset every configuration
early-stopped after 1–2 trees at ≤ +0.0002 bits validation. Its capability check
— can it recover a planted signal — is outstanding; the first two attempts were
starved by memory contention.

## 8. Status

Established:
- The feature pipeline is correct against external astronomical standards.
- The catalogue pipeline is correct and declustered.
- Model A is powerless at this dimensionality and its null means nothing.
- Under two designs, no single feature separates cases from referents.

Not established:
- Anything about the sealed test period, which remains unopened.
- Whether trees do better than a linear model, pending the capability check.
- Anything at all about slow planetary configurations: both completed designs
  are near-blind to features that barely move within a month.

**No positive result has been seen under any design.**
