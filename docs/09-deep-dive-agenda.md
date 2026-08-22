# 09 — Deep Dive Agenda: Mathematics, Equations, Methods

Target list for the next research pass. Each item states what must be derived or
obtained, the primary sources, and what it unblocks.

---

## 1. Rate-and-state response to periodic forcing ★ TOP PRIORITY

> **Second-pass update.** The simple low-pass guess is **wrong** — see
> [08-hypotheses.md](08-hypotheses.md) revision note. The literature gives a
> **two-regime** response split at T_c ≈ 2π t_a: rate tracks *stress* for T ≪ t_a,
> and *stressing rate* (amplitude ∝ 1/T) for T ≫ t_a. Heimisson & Avouac's
> analytical model is explicitly **not valid near T ≈ t_a**, which is exactly the
> region we care about.
>
> **Blocking task:** obtain full text of Ader et al. (2014) *GJI* 198, 385–413 and
> Heimisson & Avouac (2020) *GRL* 47. Both were paywalled/cert-blocked in the
> second pass. Everything downstream depends on having the actual equations.

**Derive:** amplitude response `R(ω)` and phase lag `φ(ω)` of seismicity rate to
sinusoidal Coulomb stress forcing, under rate-and-state friction.

Starting point — Dieterich (1994) rate-and-state seismicity equation:

```
dR/dt = (R / t_a) · ( τ̇ / τ̇_r  −  R )        t_a = aσ / τ̇_r
```

with R the seismicity rate, `a` the direct-effect parameter, σ effective normal
stress, τ̇_r the background stressing rate, t_a the characteristic relaxation
(nucleation) time.

For small sinusoidal perturbation `Δτ = A sin(ωt)` superposed on constant τ̇_r,
linearise and solve for the steady periodic response. Expected form: a first-order
low-pass filter with corner at ω ≈ 1/t_a, giving

```
R(ω) ∝ 1 / sqrt(1 + (ω t_a)²)          φ(ω) = −arctan(ω t_a)
```

**Verify this against the sources rather than assuming it** — the nonlinearity and
the finite-amplitude regime matter, and Beeler & Lockner identify *two* response
modes depending on stressing frequency.

**Sources:** Dieterich (1994) *JGR* 99; Beeler & Lockner (2003) *JGR* 108(B8);
Ader, Avouac et al. on tidal modulation and rate-and-state; Ader et al. (2014) on
periodic stressing of rate-and-state faults.

**Unblocks:** hypotheses 1 and 2 in [08-hypotheses.md](08-hypotheses.md) — the
spectrum inversion and the phase-vs-frequency artifact test. Highest-leverage item
on this list.

**Key questions to answer:**
- Exact form of R(ω), φ(ω), and their validity range
- How t_a maps to observable fault properties
- What the two Beeler–Lockner response modes are and where the boundary sits
- Whether phase lag is resolvable given realistic event counts

---

## 2. Tidal potential expansion and constituent catalogue

**Obtain:** the **Cartwright–Tayler–Edden** (CTE) harmonic development of the
tide-generating potential — constituent frequencies, Doodson numbers, amplitudes.
Consider also the Hartmann–Wenzel (HW95) expansion, which is more complete and
includes planetary terms directly.

**Derive:** the degree-2 (and lunar degree-3) potential from first principles:

```
V(r,ψ) = (GM/d) Σ_{n≥2} (R/d)ⁿ Pₙ(cos ψ)
```

and its decomposition into the six Doodson arguments (τ, s, h, p, N′, pₛ).

**Then:** extend with planetary mean longitudes subject to the d'Alembert
constraint Σkᵢ = 0 — the "generalised Doodson expansion" of doc 02.

**Unblocks:** hypothesis 3 (amplitude prior over harmonic orders). HW95 may
already contain the planetary terms we intend to add, which would be worth knowing
before deriving them independently.

---

## 3. Solid Earth tide: Love numbers and the response

**Obtain:** Love numbers h₂, k₂, l₂ (and frequency-dependent corrections near the
diurnal band, where the free core nutation resonance matters).

**Derive:** surface displacement and **strain/stress tensor** at a point from the
tidal potential, via the Love number formalism. This is what converts potential
into the stress tensor of doc 03.

**Sources:** IERS Conventions (2010) Chapter 7 — the authoritative and complete
treatment, including the anelastic and frequency-dependent corrections.

**Unblocks:** correct tidal tensor computation. Non-negotiable for Tier A claims.

---

## 4. Ocean tidal loading: Green's functions

**Obtain:** Farrell (1972) **load Love numbers** and the corresponding Green's
functions for displacement, tilt, strain, and gravity.

**Derive:** the loading convolution

```
L(r) = ∬_ocean  G(|r − r′|) · ρ · h(r′, t)  dA′
```

with h the ocean tide height field from a global model (FES2014, TPXO9, GOT).

**Then design:** constituent-wise precomputation. For each location and each tidal
constituent, compute the loading admittance once; runtime cost then reduces to
harmonic synthesis rather than re-convolution.

**Sources:** Farrell (1972) *Rev. Geophys.*; Agnew's SPOTL documentation;
FES2014 / TPXO model docs.

**Unblocks:** the Cochran-type high-amplitude settings where the effect is largest.
Flagged in hypothesis 8 as the likely compute bottleneck and schedule risk.

---

## 5. Coulomb stress projection

**Derive:** from the tidal stress tensor **σ** and fault geometry (strike φ_s, dip δ,
rake λ), the resolved quantities:

```
n̂ = fault normal from (φ_s, δ)
ŝ = slip direction from (φ_s, δ, λ)
Δτ  = ŝ · (σ n̂)          shear on the fault in the slip direction
Δσₙ = n̂ · (σ n̂)          normal stress
ΔCFS = Δτ + μ′ Δσₙ
```

Include the sign convention explicitly (unclamping positive) — sign errors here are
common and silently invert results.

**Also derive:** `dΔCFS/dt`, required by rate-and-state.

**Sources:** standard; King, Stein & Lin (1994) for CFS conventions. GCMT for
mechanism parameters.

---

## 6. Statistics: generalised Schuster test

**Derive:** the null distribution of the Schuster statistic at harmonic order *n*.

For N events with phases θᵢ, the order-n Schuster walk:

```
Dₙ² = ( Σᵢ cos nθᵢ )² + ( Σᵢ sin nθᵢ )²
```

Classical result: under the null, `p = exp(−Dₙ²/N)` for large N. **Confirm this
holds for arbitrary n and for non-uniformly sampled phases** — our phases come from
event times, which are not uniformly distributed in tidal phase because of
catalogue structure.

**Then:** combine across orders with Benjamini–Hochberg FDR, accounting for the
correlation between orders induced by the finite sample.

**Sources:** Schuster (1897); Tanaka et al. (2002) for seismological application;
Emter (1997) review.

**Unblocks:** analytic p-values per harmonic order, and hence a principled
multiple-testing correction — the direct antidote to the Holt & Newman failure.

---

## 7. Point process likelihoods

**Obtain and write out exactly:**

**ETAS**, space–time–magnitude form:
```
λ(t,x,y,m) = μ(x,y) f(m) + Σ_{tᵢ<t} K e^{α(mᵢ−m₀)} (t−tᵢ+c)^{−p} g(x−xᵢ, y−yᵢ; mᵢ)
```
with the standard Omori–Utsu temporal kernel and a spatial kernel g. Fix
parameterisation and estimation method (MLE vs. EM vs. Bayesian).

**Full log-likelihood** with the integral term:
```
log L = Σᵢ log λ(tᵢ,xᵢ,yᵢ,mᵢ) − ∫∫∫∫ λ(t,x,y,m) dm dx dy dt
```

**With detection function** (hypothesis 5): replace the intensity with
`λ · P_det(m,x,t)` and derive the modified integral term.

**Sources:** Ogata (1988, 1998); Zhuang et al. on stochastic declustering;
EarthquakeNPP (arXiv:2410.08226) for reference implementations.

---

## 8. Spatial machinery

- **HEALPix**: pixelisation scheme, `nside` → resolution mapping, neighbour
  finding, and the spherical harmonic transform pair (`map2alm` / `alm2map`).
- **SH synthesis** of the degree-2/3 tidal field onto the mesh — confirm the
  normalisation convention and record it (4π vs. Schmidt semi-normalised;
  mismatches are a classic silent bug).
- **Resolution kernel** derivation for the β-field inversion (hypothesis 4).

---

## 9. Numerical and precision methods

- **Clenshaw recurrence** for batched Chebyshev evaluation of SPK ephemeris segments
- **Angle-sum recurrence** for cos/sin harmonic ladders — O(N) multiply-adds
- **Two-part Julian date** arithmetic for long-baseline stability
- **Newton refinement** for angle-domain root finding (doc 06 §2): convergence
  criteria, guard against missed roots near turning points
- Error budget: what precision is actually needed? Tidal phase to ~1° is likely
  ample; ephemeris precision is almost certainly not the limiting factor, and
  confirming that would prevent over-engineering the RustSPICE layer.

---

## 10. Datasets to characterise

| Dataset | What to establish |
|---|---|
| USGS/ComCat | Full-pull mechanics, Mc history by region, format |
| GCMT | Coverage, mechanism uncertainty, magnitude threshold |
| Apollo PSE (Nakamura) | Event classification, nest locations, prior periodicity results |
| Cascadia tremor | Catalogue source, completeness, temporal coverage |
| Parkfield LFE | Availability, event counts |
| FES2014 / TPXO | Constituent coverage, resolution, licence |
| IERS EOP | Polar motion, LOD series format |
| Hydrological loading | GRACE/GRACE-FO, GLDAS — for the annual confounder |

---

## Suggested order

```
1  Rate-and-state transfer function      BLOCKED on paper access — unblock first
11 PTA red-noise statistical machinery   replaces/extends item 6; confirmed viable
12 Hydrological forcing amplitudes       hypothesis 11 — GRACE/GLDAS
6  Generalised Schuster null              baseline statistics
2  CTE / HW95 constituent catalogue       unblocks hypothesis 3
5  Coulomb projection                     needed for any real feature
3  Love numbers / IERS Ch.7               needed for correct tensors
7  Point process likelihoods              needed for any modelling
4  Ocean loading Green's functions        the expensive one
8,9,10  Machinery and data                parallel with the above
```

---

## 11. Pulsar timing array statistical machinery *(added second pass)*

**Obtain and adapt:** red noise as a Fourier-domain Gaussian process with power-law
PSD; analytic marginalisation over Fourier component amplitudes while fitting
amplitude and spectral index; hierarchical Bayesian hyperparameter marginalisation;
empirical false-alarm estimation.

**Map to our problem:** catalogue red noise (aftershock clustering, Mc drift,
hydrological loading) → fitted power-law GP. Tidal harmonics → deterministic signal
on top. Strictly more rigorous than the Schuster test, which assumes independence.

**Sources:** `enterprise` package; EPTA DR2 noise-model paper; MeerKAT PTA DR.

**Note:** PTA's empirical false-alarm estimation is the same instinct as our
time-shifted null (doc 04 §6a) — worth checking whether their machinery subsumes it.

---

## 12. Hydrological loading amplitudes *(added second pass)*

**Obtain:** GRACE / GRACE-FO gravity-derived loading, GLDAS land-surface models,
well level and snow load data. Establish amplitude, phase, and uncertainty of the
annual loading signal per region.

**Purpose:** hypothesis 11 — use hydrological forcing as a **second probe** of the
transfer function at annual period, where tidal amplitude is weak, rather than
merely subtracting it as a confounder.

**Watch:** GRACE spatial resolution (~300 km) is coarse relative to fault-scale
processes. Quantify whether that is adequate before relying on it.

Items 1, 6 and 2 are the critical path. Together they turn the harmonic framework
into a physically-grounded inverse problem with analytic significance testing —
the strategic goal identified in [08-hypotheses.md](08-hypotheses.md).
