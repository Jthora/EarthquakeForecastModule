# EarthquakeForecastModule — Handoff

Written to be read cold. Source of record for everything below is the `docs/`
tree of [PlanetaryHarmonicsModule](https://github.com/Jthora/PlanetaryHarmonicsModule),
particularly `docs/07-research-log.md`.

---

## 1. What this repository is

A probabilistic earthquake forecasting system testing whether celestial, tidal and
gravimetric features improve earthquake-rate estimates over established baselines.

**What it is not:** an earthquake *prediction* system. In seismology that term
denotes deterministic time/place/magnitude claims and is largely discredited.
*Forecasting* means probabilistic rate estimation — what CSEP evaluates. Use
"forecast" everywhere: code, docs, outputs, commit messages.

**It does not depend on AstrologyCore.** The chain forks:

```text
RustSPICE
  └─> PlanetaryHarmonicsModule
        ├─> AstrologyCore ──> Cosmic Cypher, Star Seer, Resonant Finder
        └─> EarthquakeForecastModule          ← you are here
```

Deliberate. Depending on the interpretive layer would inherit its baggage and
reviewers would be right to discount the work. **The dependency graph is part of
the argument — do not add that edge.**

---

## 2. State of the science

### Established, replicated, ours

Tidal Coulomb stress modulates slow seismicity, measured at two independent sites.

| Constituent | Period (d) | Parkfield ratio | Cascadia ratio |
|---|---|---|---|
| **M2** | 0.5175 | **223×** | **100×** |
| **N2** | 0.5274 | **45.7×** | **7.8×** |
| **O1** | 1.0758 | **137×** | **51×** |
| Mf–Sa | 13.7–365 | 0.6–2.9 (ns) | 0.4–2.7 (ns) |

Ratios are observed Schuster power over the block-shift null median; the null
expectation is 1. All three significant constituents reach the p-floor at both
sites. **8/9 constituents give the same verdict across sites.**

The sites differ in tectonic setting (strike-slip transform vs subduction
megathrust), geography (~1000 km), epoch, and **detection method** (template
matching vs envelope cross-correlation). A shared instrumental artifact would have
to survive all four differences.

**Also established:**

- **Amplitude law.** Response rises with forcing amplitude faster than linearly:
  log-log slope **3.56** against 2 for linear, surviving a re-binned null at
  p = 0.0095. Consistent with `R = R₀exp(S_T/Aσ₀)` in its non-linear regime.
- **Elastic calibration.** The M2 solid Earth tide computes to **595 Pa**, against
  Thomas et al.'s independently-inferred Parkfield `Aσ₀ = 600 Pa` — matching to 1%.
  So `S_T/Aσ₀ ≈ 1`, which independently explains the non-linear slope.
- **Response per unit stress is flat** to within ~3× from 0.5 d to 27 d. There is
  **no measured band limit** at Parkfield.

### The central untested claim

Two timescales bound the responsive band:

```text
t_n   nucleation duration   damps response above 1/t_n
      Beeler & Lockner extrapolate t_n >= 1 year for the San Andreas
T_a = 2π Aσ₀ / τ̇            critical period; response peaks here
      Ader et al. 2014 -- roughly 20-200 yr for ordinary crust

predicted responsive band for ORDINARY CRUST:   ~1 year to ~200 years
```

If true, every semidiurnal, diurnal, fortnightly and monthly constituent should be
**damped** in ordinary crust, with response appearing at Sa, the 18.61 yr nodal
term, and decadal LOD.

**Parkfield and Cascadia are the control, not the test.** Both are *tremor*, where
`T_a` is short — so short-period response there is *expected* and says nothing
about ordinary crust.

> **P3.4 is the whole point of this repository.** Measure `R(ω)` for ordinary
> crust and compare. The same shape as tremor refutes the band prediction; the
> mirror image confirms it. Everything else here is instrument.

⚠ `t_n ≥ 1 yr` is a lab-to-field extrapolation — the paper's own inference, not a
measurement. Carry it with wide uncertainty.

### Sample size

Beeler & Lockner 2003 (*JGR* 108(B8), **free from USGS**), equation 18:

```text
N >= ln(P_rw) / (Δτ_u / (2 a σ_n))²
```

Worked examples give **6.2×10³–5.5×10⁴** events; the abstract states daily Earth
tides need **">13,000 earthquakes to detect."** N scales as the **inverse square**
of normalised stress amplitude, so 10× amplitude cuts the requirement 100×.

*(An early draft of our notes said 10⁵–10⁶, taken from a search summary and never
verified. Wrong by two orders of magnitude.)*

### The invariant that constrains every claim

**⟨R⟩ = r exactly** (Heimisson & Avouac 2020, eq. 6). Oscillatory stress does not
change the mean rate — only the timing. **Tides redistribute *when* events occur;
they do not create them.** Any model, code path or output implying otherwise is a
bug.

---

## 3. What PlanetaryHarmonics provides

Rust crate `ph-core`, **72 tests**, in the submodule at `modules/`. Consume it; do
not reimplement.

| Module | Provides |
|---|---|
| `tidal` | Degree-2 tensor `T_ij = (GM/d³)(3n̂ᵢn̂ⱼ − δᵢⱼ)`, superposition, eigendecomposition, principal axis |
| `field` | Tidal fields from real ephemeris — Earth in ITRF/IAU_EARTH, Moon in MOON_PA |
| `fault` | Aki & Richards geometry, traction decomposition, `ΔCFS = τ + μ′σₙ`, rotation to local NED, **linear Coulomb coefficients** |
| `love` | Elastic response: tensors (s⁻²) → stress (Pa); `T_a` relations |
| `doodson` | **Analytic constituent phases** from the six fundamental arguments |
| `stats` | Generalised Schuster, periodogram, time-shifted and block-shift nulls |
| `phase` | Tidal phase from a sampled quasi-periodic forcing |
| `demod` | Complex demodulation (see trap 5 before using its phase output) |
| `catalog`, `apollo`, `parkfield`, `cascadia` | Event types and ingestion |

**Validated end to end.** Phase 1 recovered five known deep moonquake
periodicities from the Apollo catalogue to better than 0.21%.

### Two consumption paths

**CLI** — `ph-features`, zero dependencies, emits CSV with a full provenance
header (frame, epoch system, aberration, kernels, site, geometry, elastic
constants, tiers, units, caveats):

```bash
ph-features --lat 35.635 --lon -120.150 \
            --strike 137 --dip 90 --rake 180 --mu 0.4 \
            --start 2001-01-01 --days 8400 --step 0.02 \
            --out features.csv
```

**Python** — `maturin build --release` in `crates/ph-py`, then install the wheel.
Exposes constituent phases and periods, stress/strain conversion, Schuster, and
the block-shift null. Verified on Python 3.14 via a cp39-abi3 wheel.

Ephemeris and tensor computation are deliberately **not** in the Python bindings —
a SPICE session is not `Send`, and that work is batch. Heavy geometry through the
CLI, analysis primitives through the bindings.

**Do not strip the provenance header.** A feature whose frame or epoch system is
ambiguous is worthless in a statistical test, and across a repo boundary that
ambiguity is easy to introduce.

---

## 4. What this repository must build

1. **USGS ComCat ingestion** + magnitude-of-completeness per region and epoch
2. **GCMT focal mechanisms** → fault geometry. Note this is the *easy* case: on
   Earth the mechanism is known, so no orientation search is needed (the Moon
   required one, and `ph-core::fault` supports both)
3. **ETAS baseline**, fitted and **frozen**
4. **Residual model** — `λ = λ_ETAS · exp(f_θ(features))`
5. **`R(ω)` for ordinary crust** ← **P3.4, the point**
6. **β(x,t) sensitivity field**
7. **CSEP evaluation** via `pyCSEP`

---

## 5. Methodology — non-negotiable

**It is a point process, not classification.** Target the conditional intensity
`λ(x,y,t,m)` on point-process log-likelihood. The integral term makes the *absence*
of earthquakes informative. Classification gives 99.9% accuracy and learns nothing.

**Fit ETAS first, freeze it, learn only the residual.** With a prior pulling
`f_θ → 0`, the model becomes structurally incapable of taking credit for clustering
ETAS already explains, and `exp(f_θ)` is directly interpretable.

**Mc filtering is mandatory.** Detection improved dramatically over any long
catalogue; uncorrected, that trend projects onto long-period features and
manufactures signal — fatal given that the band prediction points at exactly those
periods.

**FDR on any feature scan.** Benjamini-Hochberg, or pre-register a small set.

**Report information gain per earthquake** in bits/event via the CSEP T-test. Run
N/S/M/L tests. Do not invent metrics.

### The standing rule, earned six times

> **Specify the null before running, and state what structure it preserves.**

Every trap below came from choosing a null after seeing the data. The first test
run with the null fixed in advance was also the first result to survive
multiple-comparison correction. Not a coincidence.

---

## 6. The traps — read this section twice

All six share one shape: **a statistic that looks decisive while silently answering
a different question.**

**1 — Time-shift null degenerate when the catalogue shares the forcing's period.**
Deep moonquakes are locked near the anomalistic month; so is the forcing. A global
shift *rotates* the phase cluster without diluting it, so `D²` is near-invariant.
Analytic Schuster gave p = 10⁻⁸⁹; the empirical null gave **p = 0.70**. An 88-order
overstatement.

**2 — Uniform-time null tests temporal clustering, not tidal alignment.** Drawing
random event times destroys the catalogue's clustering, so the test silently
becomes "are events clustered in time at all?" — trivially yes. Gave **73/74**
moonquake nests significant; a shift-preserving null gave **17/74**, and **0/74**
after FDR.

**3 — Raw period folding measures the detector.** On a detection-limited catalogue,
folding event times on trial periods finds the instrument. At Parkfield, S1
(exactly 24.000 h, essentially no body-tide amplitude) reaches Schuster power
**16,245** against an expectation of 1. **K1 (23.93 h) sits at 1.16× that floor and
is unusable.** S2 (12.000 h) exceeds M2 by 4.2×, which is backwards for a body tide.

**4 — Per-bin nulls compromised when the binning variable derives from the
forcing.** Binning by tidal amplitude inherits spring–neap structure. The
highest-amplitude bin had the *highest* response and yet p = 0.45. The fix is to
null the *claim* (the trend), by re-binning at shifted times.

**5 — Time-shift null degenerate against a single demodulated constituent.** Trap 1
recurring. Demodulation makes a band a near-pure tone, against which a shift is a
pure rotation. **This was documented in the module and walked into anyway.**
Documenting a trap does not confer immunity; it reappears wearing different
clothes.

**6 — Sham-frequency nulls return the leaking constituent's phase.** Running the
procedure at tide-free frequencies gave `D²/N ≈ N` — every event at one phase,
where there is no tide. Where there is no genuine power, `z̄` is dominated by
leakage at ω′, so the reported phase collapses to `ω′t`. A tide-free frequency is
not a neutral baseline; it is a relabelled copy of the dominant band.

**And one that is not a null at all:** a raw response measurement is not a transfer
function. Dividing by the forcing amplitude is not a refinement, it *is* the
definition. Skipping it produced a confident, wrong, physically-flavoured
conclusion ("the response is band-limited") that survived a full write-up before
normalisation caught it. Every statistic involved was sound.

### The null that works for a single constituent

A **global** shift can never work — `D²` is rotation-invariant and a global shift
*is* a rotation. Shift **each block independently**, block length
`max(4 × period, 30 d)`: within-block clustering preserved, between-block alignment
randomised. Implemented in `ph-core::stats` and exposed to Python.

---

## 7. Validity gates — run before believing any positive result

Cheap, need no external data, and they are gates rather than analyses.

**M2 vs S2.** S2 is exactly 12.000 h, locked to the day–night cycle, and carries
the solar *thermal* tide. M2 is 12.42 h and precesses through local solar time,
decorrelating from time-of-day artifacts. Signal at S2 but not M2 → artifact.
Likewise **K1 (23.93 h) and S1 (24.00 h) are unusable**; O1 (25.82 h) and P1
(24.07 h) are safe.

**Alias analysis.** Enumerate catalogue periodicities (daily, weekly, seasonal
maintenance, network upgrades), compute beats against the constituent list,
blacklist collisions.

**Two free amplitude knobs.** Lunar distance varies 5.5% over the anomalistic
month and tides go as 1/d³ → **18% amplitude modulation** at 27.55 d. The 18.61 yr
nodal cycle modulates diurnal amplitudes ~±11%, giving **seven envelope cycles**
over a 130-year catalogue.

**Same-band constituent ratios.** Within a band the transfer function is roughly
constant, so response ratios should follow amplitude². At Parkfield **O1/Q1 = 31.6
against 28 predicted** (good) while **M2/N2 = 11.5 against 28** (off by 2.4×) —
unexplained, and worth understanding rather than glossing.

**Hydrological loading is a confounder *and* an instrument.** Seasonal groundwater,
snow and atmospheric loading are annual, and so are Sa and Ssa. Any annual
"celestial" correlation is confounded by default. But hydrological amplitude is
independently measurable from GRACE/GRACE-FO and GLDAS, so model it jointly: it
becomes a **second probe of the transfer function** at annual period, where tidal
amplitude is weak. Precedent: *Science Advances* sciadv.ady6350.

---

## 8. Data — all free, no credentials

| Source | Contents |
|---|---|
| USGS ComCat | Global earthquake catalogue, public API |
| GCMT | Focal mechanisms |
| IRIS / EarthScope | Waveforms and catalogues |
| IERS EOP | Polar motion, length of day |
| GRACE / GRACE-FO, GLDAS | Hydrological loading |
| FES2014 / TPXO9 / GOT | Ocean tide models |
| NAIF | DE440/DE441 SPICE kernels |

**Watch for silent truncation.** The PNSN tremor API caps responses at 20,000
events with HTTP 200 and no truncation flag; a yearly request returned exactly
20,000 and produced a tidy-looking catalogue that was quietly wrong. Assume any
API paginates or caps until proven otherwise, and check for suspiciously round
counts.

**Literature without institutional access:** author self-archived pages first
(highest yield), then **USGS Publications Warehouse** — USGS-authored work is public
domain and covers much of this field (Beeler, Cochran, Hardebeck) — then arXiv,
**ESS Open Archive** (AGU preprints: GRL, JGR), EarthArXiv, and the Unpaywall API
by DOI.

---

## 9. Prior art — use it, do not rebuild

- **RECAST** — neural temporal point process; matches or beats ETAS on Southern
  California given enough events. `github.com/keliankaz/recast`
- **EarthquakeNPP** — benchmark suite, arXiv:2410.08226
- **pyCSEP** — the evaluation framework, `cseptesting.org`
- **Pulsar timing array methods** — red noise as a Fourier-domain Gaussian process,
  analytic marginalisation over Fourier amplitudes, empirical false-alarm
  estimation. Same problem shape (small periodic signal, strongly red noise) and
  far more rigorous than the Schuster test, which assumes independence. Tooling:
  `enterprise`.

---

## 10. Terminology

| Avoid | Use |
|---|---|
| prediction | probabilistic forecast, conditional rate estimate |
| nexus event | coherence maximum, commensurability alignment |
| resonance point | phase coherence peak |
| astrology feature | derived planetary feature |
| base-N harmonics | harmonic order N of relative longitude |

Label every feature by tier: **A** established physics, **B** established method /
novel application, **C** exploratory. Report Tier C separately, never as mechanism.

---

## 11. Open questions

- **Does the 1 yr – 200 yr band prediction survive measurement for ordinary
  crust?** The project's central question, not answerable from the literature.
- **M2/N2 responds 2.4× more than the amplitude law predicts.** Real frequency
  dependence within the semidiurnal band, or an unfound leak?
- **`T_a` is not located anywhere.** Parkfield's `R(ω)` is flat across the
  measurable range, so no peak has been seen at either site.
- **Long-period power.** Mf through Sa are non-significant at both tremor sites, but
  that is a power limitation, not evidence of absence. Only a long ordinary-crust
  catalogue can resolve it.
- **What publishes if the answer is null?** **Decide before running P3.4.** The
  honest answer is that a measured transfer function for ordinary crust is
  publishable *whatever shape it has*, alongside the moonquake and two-site tremor
  validations. Pre-committing removes the incentive to keep slicing until something
  crosses p < 0.05.
