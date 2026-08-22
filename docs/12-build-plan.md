# 12 — Build Plan

**Decision (2026-08-21): start building.** The central open question — whether the
responsive band is ~1 yr to ~200 yr (doc 07, fourth pass) — is not resolvable by
further reading. The literature contains no global, multi-decade response spectrum;
that gap *is* our novelty. It requires measurement.

Four research passes each corrected the previous, but the corrections moved from
structural (pass 2: wrong model shape) to parametric (pass 4: wrong constant). That
is diminishing returns.

---

## Two things that unblocked without announcing themselves

**1. HW95 is not a blocker.** The catalogue is a *frequency-domain* decomposition.
We can compute the tidal potential and tensor **directly in the time domain from
DE440 positions**, then do our own harmonic decomposition by least squares on the
series we generated. HW95 is a cross-check and an amplitude prior — valuable for
hypothesis 13a's ratio tests, but it gates nothing.

**2. The band prediction defers the most expensive component.** Ocean tidal loading
— the Farrell Green's-function convolution flagged as the schedule risk (doc 08 §8)
— is dominated by semidiurnal and diurnal constituents. If the responsive band is
years-to-decades, its contribution there is small relative to solid tide. The
prediction that reshapes the science also removes the worst engineering from the
critical path.

**And the measurement does not depend on the unsettled theory.** Heimisson &
Avouac's model being invalid near T ≈ t_a does not matter: measuring a response
spectrum is model-free. Theory interprets the spectrum; it does not produce it.

---

## Phase 1 — Deep moonquake validation

**Not the earthquake system.** Phase 1 is simultaneously the known-answer test and
the minimal end-to-end pipeline.

Why moonquakes first:

- ~10⁴ events — right-sized
- Tidal forcing is **dominant**, not marginal — there is a known answer to check
- No ocean loading, no hydrology, no Mc drift, no cultural noise
- Exercises the whole chain: ephemeris → tensor → phase → spectrum → statistics

**Success criterion:** recover the known deep moonquake periodicities — monthly
(27.55 d), 206-day, and ~6-year — from the Apollo PSE catalogue, with a
time-shifted null showing they are not artifacts.

**If it fails, everything downstream is worthless and we learn it cheaply.**

### Scope

| # | Component | Notes |
|---|---|---|
| 1 | Workspace + `ph-core` crate consuming `rustspice-core` | Rust library, no WASM |
| 2 | Ephemeris layer — batched geometric states | `Aberration::None` (doc 10 §2) |
| 3 | Tidal tensor — degree 2, plus lunar degree 3 | `T_ij = (GM/d³)(3n̂ᵢn̂ⱼ − δᵢⱼ)` |
| 4 | Principal axis / concentration node | Symmetric 3×3 eigendecomposition |
| 5 | Harmonic encoding + time-domain decomposition | No catalogue needed |
| 6 | Generalised Schuster spectrum at arbitrary order | `p = exp(−D²ₙ/N)` |
| 7 | Time-shift null distribution | Doc 04 §6a |
| 8 | Apollo PSE ingestion | NASA PDS, public |
| 9 | Phase 1 validation run | Against known periodicities |

Every component is reused by Phases 2 and 3.

### Kernels for Phase 1

`naif0012.tls` (leap seconds), `de440.bsp` (ephemeris), `pck00011.tpc` (constants),
**`gm_de440.tpc` (GM values — the tensor is GM/d³)**, and `moon_pa_de440_*.bpc` +
`moon_de440_*.tf` for selenographic frames.

---

## Phase 2 — Tectonic tremor

High-SNR terrestrial testbed. Tremor's T_a sits inside the tidal band (doc 07,
third pass), so response is strong and clearly frequency-dependent.

Adds: Cascadia / Parkfield LFE catalogue ingestion; the **full exponential response
form** `R = R₀exp(S_T/Aσ₀)` with M > 1 carried, since Aσ₀ ≈ 6×10⁻⁴ MPa puts tremor
firmly in the non-linear regime.

## Phase 3 — Earthquakes

Adds: USGS ComCat ingestion, Mc analysis, GCMT focal mechanisms, Coulomb projection
onto fault geometry, ETAS baseline, the β(x,t) field, hydrological covariates, CSEP
evaluation.

**Validity gates from hypothesis 13 run before any positive result is believed:**
M2-vs-S2 cultural-noise discrimination, explicit alias analysis, and the
perigee–apogee and nodal-envelope amplitude tests.

---

## Explicitly deferred

| Deferred | Why |
|---|---|
| Ocean loading convolution | Band prediction says it matters little in the responsive range. Revisit if short-period response appears. |
| Full ETAS-residual ML stack | Needs the spectrum measurement first, to know which features matter. |
| Multi-basis feature explosion | Let the HW95 amplitude prior constrain it rather than generating everything. |
| WASM / TypeScript surface | For Cosmic Cypher, Star Seer, Resonant Finder. Science is not settled enough to freeze an API, and Rust-library composition means waiting costs nothing. |
| Spherical harmonic / HEALPix global synthesis | Phase 3. Phase 1 is point-wise. |

---

## Finish before the Phase 3 split

Two physics pieces belong in `ph-core`, not in EarthquakeForecastModule, and
building them here means the forecasting repo starts with a complete feature layer
rather than a half-built one:

1. **Ephemeris-based tidal phase.** Compute phase from real Earth–Moon–Sun geometry
   rather than folding on trial periods. Required by the time-shift null, which is
   degenerate against a single exact frequency (see `ph-core::stats`). Needs the
   lunar kernels.
2. **Coulomb projection.** `ΔCFS = Δτ + μ′Δσₙ` and `dΔCFS/dt` from the tidal tensor
   and fault geometry (strike/dip/rake). This is physics, it is shared, and
   forecasting needs it on day one.

Also plan **`pyo3`/`maturin` bindings** — needed by Layer 2's own transfer-function
ML regardless of any split (`docs/13-ml-stack.md` §2c), and the natural interface
for a Python forecasting stack.

## Riding-along research task

Doc 09 §10 still lists the Apollo PSE catalogue as "to characterise": event
classification scheme, nest locations, and the published periodicity results we are
validating against. An afternoon, not a research phase — do it alongside the build.

Also still open, but non-blocking: free access to HW95/ETERNA catalogue data
(doc 07, fourth pass). The Kudryavtsev & Cionco (2025) release ships in HW95 format
and is the most promising free route.

---

## Guiding constraint

From doc 00's claim discipline, reinforced by Heimisson & Avouac eq. 6:

> **⟨R⟩ = r exactly.** Oscillatory stress does not change the mean rate — only the
> timing. Tides redistribute *when* events occur; they do not create them.

Any code path or output that implies otherwise is a bug.
