# 00 — Project Framing and Claim Discipline

## What this module is

A computational framework for deriving **astronomical, tidal, gravimetric, and
multi-basis harmonic features** from high-precision ephemerides, for statistical
analysis of Earth-system events.

Formal descriptor:

> Celestial–Terrestrial Harmonic Feature Analysis for Probabilistic Earthquake Forecasting

## What this module is not

It does not interpret. It does not assign meaning. It emits numbers with
documented units and provenance. Every output is a deterministic function of
time, position, and body ephemerides.

The scientific question is **"does this feature set improve forecast skill over an
established baseline?"** — not "do the planets cause earthquakes."

## Claim discipline

Three tiers, always labelled distinctly in code, docs, and outputs:

**Tier A — Established physics.** Solid Earth tides, ocean tidal loading, Coulomb
failure stress, pole tide, rate-and-state friction. Peer-reviewed, quantified,
uncontested mechanism. Tidal triggering of earthquakes is *real and small*
(see [01-literature.md](01-literature.md)).

**Tier B — Established statistics, novel application.** Generalised Doodson
expansion including planetary arguments; Fourier harmonic encodings of relative
angle; multi-body commensurability terms. The mathematics is standard celestial
mechanics; applying it as an ML feature basis is the new part.

**Tier C — Exploratory.** Anything without a proposed physical mechanism.
Permitted, but must be pre-registered, FDR-corrected, and reported separately
from Tier A/B results. Never presented as mechanism.

## Terminology

Internal vocabulary is generative and can stay. External-facing work uses:

| Internal | External |
|---|---|
| prediction | probabilistic forecast / conditional rate estimate |
| nexus event | coherence maximum / commensurability alignment |
| resonance point | phase coherence peak |
| timestream | temporal grid / evaluation epoch series |
| astrology feature | derived planetary feature |
| base-N harmonics | harmonic order N of relative longitude |

Avoid calling every base-number representation a "harmonic." Reserve:
- **Fourier harmonics** — periodic frequency components
- **Spherical harmonics** — global spatial fields
- **Tidal constituents** — astronomical tide frequencies (Doodson-indexed)
- **Orbital resonance terms** — integer frequency commensurabilities
- **Modular / radix encodings** — the generalised base-number representation

## The cautionary example

Holt & Newman (2025), *Tidal Triggering of Magnitude 7+ Earthquakes by Jupiter*
(arXiv:2508.07064), ran **28,782 chi-squared tests** over interval lengths of
5–45 days and reported **1,071 statistically significant** intervals.

At α = 0.05, chance alone predicts ≈ 1,439 significant results. They found
**fewer significant results than the null hypothesis expects.**

(The tests are heavily correlated — overlapping intervals — so the effective
number of independent tests is well below 28,782 and the naive expectation is
not exactly 1,439. The conclusion stands regardless: without multiple-comparison
correction and without accounting for test correlation, the reported count is
not distinguishable from noise.)

This is the precise failure mode this project must be built to avoid. It is also
the reason for the two hard requirements in
[04-ml-architecture.md](04-ml-architecture.md): FDR correction on any feature
scan, and a **time-shifted null distribution** rather than random shuffling.
