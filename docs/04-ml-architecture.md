# 04 — ML Architecture and Validation

## 1. This is a point process, not a classification problem

The failure mode that kills celestial-correlation studies is framing the task as
"label each grid cell × time bin as quake / no-quake." You get 99.9% accuracy,
learn nothing, and class imbalance hides everything.

The correct target is a **conditional intensity function**

```
λ(x, y, t, m)  —  expected event rate per unit space × time × magnitude
```

trained by maximising **point-process log-likelihood**:

```
log L = Σ_events log λ(xᵢ, tᵢ, mᵢ)  −  ∫∫∫ λ(x, t, m) dx dt dm
```

The integral term makes the *absence* of earthquakes informative. That is where
most of the information lives, and the classification framing discards it.

## 2. Baseline: ETAS, non-negotiable

The baseline must be **ETAS** (Epidemic-Type Aftershock Sequence) — the operational
standard and the reference CSEP tests against. Plus a time-independent background
from smoothed seismicity.

If harmonic features do not beat ETAS, there is no result. If they do, it is
publishable regardless of mechanism.

## 3. Architecture: learn the residual, never train alongside

```
λ(x, t, m) = λ_ETAS(x, t, m) · exp( f_θ(harmonic features) )
```

**Fit ETAS first. Freeze it.** Then learn only the multiplicative modulation
`f_θ`, with a prior pulling it toward zero.

Why this is the key design decision:

1. The model is **structurally incapable** of taking credit for clustering ETAS
   already explains.
2. `exp(f_θ)` is **directly interpretable**: "planetary features modulate rate by
   ±2.8%" is a real, quotable, falsifiable claim.
3. It matches the physics — tidal triggering *is* a small multiplicative rate
   modulation on a background process. The architecture encodes the hypothesis.

This single choice does more for credibility than any validation bolted on later.

## 4. The stress-state target (from Beaucé / Ide)

Per doc 01, the stronger research target is not rate modulation but
**time-varying tidal sensitivity** as a proxy for proximity to failure.

Parameterise the coupling coefficient itself as a learned field:

```
λ = λ_ETAS · exp( β(x, t) · ΔCFS_tidal(x, t) )
```

where **β(x, t)** — the sensitivity — is the learned output, smoothed in space and
time. Then:

- β rising in a region is the Beaucé precursor signal, generalised
- β is a physically interpretable scalar, not a black-box score
- The forecast is the model; β is the diagnostic

This is the most defensible and most novel modelling contribution available.
Recommend making it the primary target.

Secondary target from Ide et al.: condition the **magnitude distribution** (b-value)
on tidal amplitude, not just the rate.

## 5. Model tiers

| Tier | Model | Purpose |
|---|---|---|
| 0 | ETAS + smoothed background | Baseline. Must beat this. |
| 1 | ETAS × exp(GLM on Fourier features) | Interpretable; coefficients = Schuster stats |
| 2 | ETAS × exp(NN residual) | Nonlinear interactions; still bounded |
| 3 | β(x,t) sensitivity field | Stress-state proxy; primary research target |
| 4 | Neural temporal point process | Full flexibility |

**Tier 4 prior art:** RECAST (Kaz et al.) — encoder–decoder GRU neural temporal
point process, matches or beats ETAS on Southern California given a large enough
catalogue. Code: `github.com/keliankaz/recast`. Benchmark: **EarthquakeNPP**
(arXiv:2410.08226). Use these rather than building from scratch.

**Spatial:** use **HEALPix** for the global mesh — equal-area and natively
compatible with spherical-harmonic synthesis, which is how the tidal field is
represented anyway (doc 03). H3 is the alternative but lacks the SH affinity.

## 6. Validation — the part that determines whether this is science

### 6a. Time-shifted null distribution

**Do not randomly shuffle.** Shuffling destroys aftershock autocorrelation and
gives falsely tight nulls. Celestial features are smooth and quasi-periodic,
making them maximally vulnerable to this.

Instead: **shift the entire catalogue in time** by large offsets (±1 to ±50 years)
and recompute the full pipeline. This preserves the catalogue's internal
clustering *and* the celestial signal's structure, and breaks **only their
alignment**. Build the null from hundreds of shifts.

Cheap, brutal, and very hard to argue with. This is the single most important
statistical safeguard in the project.

Avoid shift offsets near integer multiples of dominant periods (1 yr, 18.61 yr) —
those partially preserve alignment.

### 6b. Magnitude of completeness

Filter to events above the **region- and epoch-specific Mc**. Non-negotiable.
Detection capability improved dramatically over the catalogue's span; uncorrected,
that trend projects onto any long-period feature and manufactures signal.

### 6c. Multiple-comparison correction

Any feature scan gets **Benjamini–Hochberg FDR**. Better: pre-register a small
feature set and confirm on held-out data.

See the Holt & Newman analysis in doc 00 for what happens without this.

### 6d. Declustering

Either explicitly model aftershocks (ETAS, preferred) or decluster. Never treat
raw catalogue events as independent samples — one mainshock generates thousands of
correlated "events," and any slowly-varying feature will appear correlated with
the aftershock sequence.

### 6e. CSEP-style evaluation

Report **information gain per earthquake (bits/event)** over baseline, via the
CSEP **T-test** of equal predictive ability. Run **N-test** (number),
**S-test** (spatial), **M-test** (magnitude), **L/PL-test** (likelihood).

Use `pyCSEP` — cseptesting.org. Do not invent metrics; the field has standard ones
and using them removes an entire class of objection.

### 6f. Prospective testing

Reserve the final years. Better: register a genuine forward forecast with CSEP.
Anything purely retrospective will be doubted, correctly.

## 7. Power analysis — do this first

Beeler & Lockner: **10⁵–10⁶ events** needed for robust tidal correlation, with a
10× stress amplitude increase cutting the requirement 100×.

Before training anything, compute the **minimum detectable effect size** given
catalogue size per region. This tells you which regions are worth analysing at
all, and prevents burning months on underpowered subsets.

Expected outcome: only high-tidal-amplitude settings (subduction zones, coastal,
mid-ocean ridges) will have adequate power. That is a finding, not a limitation —
and it matches the physical prediction.

## 8. Pipeline order

1. Full USGS catalogue pull; Mc analysis per region/epoch
2. Power analysis → region selection
3. GCMT focal mechanisms → fault geometry
4. Feature generation (docs 02, 03) with provenance metadata
5. ETAS fit, frozen
6. Tier 1 GLM → Schuster-comparable coefficients
7. Time-shifted null (hundreds of shifts)
8. FDR correction
9. Tier 3 sensitivity field
10. CSEP evaluation vs. ETAS
11. Held-out / prospective confirmation
