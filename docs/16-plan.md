# 16 — Working Plan

Established 2026-08-21, after Phase 1's periodicity validation passed
([07-research-log.md](07-research-log.md)).

---

## Investigation findings that shaped this plan

**1. All Phase 1 kernels are available and small.** Confirmed at NAIF:

| Kernel | Purpose |
|---|---|
| `naif0012.tls` | Leap seconds |
| `de440s.bsp` | Ephemeris — **32.7 MB**, covers 1849–2150 (full `de440` is 114 MB and unnecessary) |
| `pck00011.tpc` | Body constants |
| `gm_de440.tpc` | GM values — the tensor is GM/d³ |
| `moon_pa_de440_200625.bpc` | Lunar orientation with libration |
| `moon_de440_250416.tf` | Lunar frame definitions |

Fetch with `./scripts/fetch-kernels.sh`. **Task 1 is unblocked.**

**2. Deep moonquake focal mechanisms are poorly constrained — and that is useful.**

Weber, Bills & Johnson (2009), *JGR* 114, E05001, and Weber & Knapmeyer
(NTRS 20120013649, **free**): Apollo's limited source/receiver geometry
*prohibits* determining fault parameters from first-motion polarities. Their
S/P-amplitude inversion reduces the fault-plane parameter space by only about half
on average; the best-constrained cluster eliminates 72% of the focal sphere.

Their method: **grid search over strike/dip/slip, with the failure criterion that a
linear combination of shear and normal stress best approximates a constant** at
event times. Data: 106 events across 25 clusters, 37 of them from A001.

**This reorders the plan.** Three consequences:

- **Coulomb projection must treat fault orientation as a *parameter*, not an
  input.** On the Moon we search for it; on Earth we fix it from GCMT. **One code
  path, two uses** — and building it against the harder case first means the Earth
  case is a simplification rather than a rewrite.
- **Weber 2009 is a published benchmark**, giving Phase 1 a *second* known-answer
  validation — this time of the stress projection, not just the timing analysis.
- They find **stress *rates*** matter for some clusters, independently supporting
  the `dΔCFS/dt` emphasis in [03-tidal-tensor.md](03-tidal-tensor.md).

So **task 4 (Coulomb projection) moves up**, and merges with task 3.

---

## Track A — Finish Phase 1  ◀ start here

| # | Task | Depends on |
|---|---|---|
| **A1** | Fetch kernels; wire `ph-core` to load them | — |
| **A2** | Compute lunar tidal tensor from real Earth + Sun geometry, in `MOON_PA` | A1 |
| **A3** | Tidal phase from the full quasi-periodic forcing (not trial-period folding) | A2 |
| ~~A4~~ | ~~Time-shift null on the pooled catalogue~~ — **does not work; see below** | — |
| **A4′** | **Per-nest** phase test: is the preferred phase consistent within a nest and structured across nests? | A3 |
| **A5** | Coulomb projection with fault orientation as a searchable parameter | A2 |
| **A6** | Validate A5 against Weber 2009 per-cluster fits; use `nakamura_2005_dm_locations.csv` (already downloaded, unused) | A5 |

**A3 is the linchpin.** The time-shift null is degenerate against a single exact
frequency, so A4 cannot be done correctly without it.

> **A4 failed as specified, 2026-08-21.** The time-shift null is degenerate
> whenever the **catalogue** shares the forcing's period — not merely when the
> forcing is a single sinusoid. Deep moonquakes are locked near the anomalistic
> month and so is the forcing, so a shift rotates the phase cluster without
> diluting it. Measured: analytic Schuster p = 1.07e-89, empirical p = 0.699.
> See [07-research-log.md](07-research-log.md).
>
> **Pooled-catalogue tests cannot produce a falsifiable Phase 1 claim.** A5/A6
> (per-nest Coulomb projection against Weber's constraints) are now the critical
> path, not a bonus validation.

**Exit criterion:** the Coulomb projection reproduces Weber's per-cluster fault
constraints, and preferred phase is consistent within nests and structured across
them.

## Track B — Prepare the split

| # | Task | Depends on |
|---|---|---|
| ~~B1~~ | **Done.** cp39-abi3 wheel, verified on Python 3.14 | ✓ |
| **B2** | CLI emitting feature vectors as CSV/Parquet with full provenance | — |
| ~~B3~~ | **Done.** [EarthquakeForecastModule](https://github.com/Jthora/EarthquakeForecastModule) created with an expanded `HANDOFF.md`, README and this repo as a submodule. Fresh recursive clone verified: full chain resolves and 72 tests pass | ✓ |

B1 is needed by Layer 2's own ML regardless of any split
([13-ml-stack.md](13-ml-stack.md) §2c), so it is not speculative.

## Track C — Phase 2, tremor

| # | Task | Depends on |
|---|---|---|
| ~~C1~~ | **Done.** Shelly Parkfield LFE catalogue: 1,528,117 events, 88 families, 23.1 yr, USGS public domain | ✓ |
| **D1′** | **Done.** M2/S2 gate fires — S2 is 4.2× M2 and K1 sits at 1.16× the S1 artifact floor. Raw period folding on this catalogue measures the detector | ✓ |
| ~~C2~~ | **Done.** ΔCFS on SAF geometry, whole-day-shift null fixed in advance. **9/12 families survive BH-FDR**; phases coherent within 71°. First result to survive correction | ✓ |
| **C2c** | Full exponential `R₀exp(S_T/Aσ₀)` with M > 1 — amplitude modelling, not timing | C2 |
| **C2b** | Alias analysis — blacklist beats against the 24 h detection cycle | C1 |
| ~~C3~~ | **Done.** Amplitude law: `D²/N` rises 71× over a 3× amplitude range, log-log slope **3.56** vs 2 for linear. Trend survives its own re-binned null (p = 0.005, null median 0.66). Steeper than linear → non-linear exponential regime | ✓ |

## Track D — Test the band prediction

| # | Task | Depends on |
|---|---|---|
| **D1** | M2-vs-S2 and alias-analysis validity gates | A3 |
| ~~D2~~ | **Closed as a blocker.** Constituent *phase* needs no catalogue — it is an analytic combination of six astronomical arguments, now in `ph_core::doodson`. Only *amplitudes* need HW95/KSM03, which F2 will want | ✓ |
| **D3** | Broadband response spectrum — tides, hydrological, pole tide, LOD in one fit | C3, D2 |

**D3 is the project's central question.** Everything before it is instrument-building.

## Track E — Layer 2 ML

| # | Task | Depends on |
|---|---|---|
| **E1** | Transfer function estimator — parametric fit first, GP second | C3, B1 |
| **E2** | Symbolic regression on the `T ≈ t_a` gap where the published theory fails | E1 |
| **E3** | Ocean loading surrogate | **only if** the band prediction fails |

E3 is explicitly conditional. If the responsive band is years-to-decades, ocean
loading barely matters and this is wasted effort.

## Track F — Research, parallel

| # | Task |
|---|---|
| **F1** | PTA red-noise machinery — port `enterprise`-style methods |
| **F2** | IERS Conventions Ch. 7 — correct Love numbers, frequency-dependent corrections |
| **F3** | Obtain Weber 2009 full text (ResearchGate / AGU PDF) for exact cluster parameters |

---

## Critical path

```
A1 → A2 → A3 → A4        Phase 1 complete
          A2 → A5 → A6   Coulomb validated against Weber
                    B1, B2 → B3        split
          A3 → C1 → C2 → C3 → D3       the central question
```

**A3 gates almost everything.** It unblocks the null (A4), tremor (C1), and the
validity gates (D1).

## Decision points

1. **After A4** — if the known periods do *not* survive the null, stop and fix the
   pipeline. Nothing downstream is meaningful.
2. **After A6** — if Coulomb projection cannot reproduce Weber's constraints, the
   stress projection is wrong and Earth work must wait.
3. **After D3** — if the band prediction fails and short periods dominate, ocean
   loading (E3) returns to the critical path and the schedule changes materially.
4. **Before any Phase 3 work** — pre-commit to what gets published if the answer is
   null: a rigorous upper bound plus the moonquake and tremor validations. Deciding
   now removes the incentive to keep scanning until something crosses p < 0.05.

---

# Revision 2 — 2026-08-22, after Tracks A and C

Tracks A and C are complete. This section supersedes the track tables above.

## Where the project stands

**Two independent positive results at Parkfield**, on different statistics with
different nulls: C2's phase clustering (9/12 families surviving Benjamini-Hochberg)
and C3's amplitude law (slope 3.56, trend survives its own re-binned null).

**Four documented traps**, each caught only by checking what the null was actually
testing. That record is now a deliverable in its own right.

**The band prediction — the project's central question — is untested**, and two
things block it.

## The two blockers, stated plainly

### 1. Everything sits at one frequency

Every result so far is at **M2, 12.42 h**. A transfer function is a *curve*;
we have measured one point on it. Worse, M2 sits deep inside the range doc 07
predicts should be **damped** for ordinary crust — Parkfield works precisely
because tremor's `T_a` is anomalously short.

We cannot say anything about the 1 yr–200 yr band from a 12.42 h measurement.

### 2. No Love numbers, so `T_a` is unlocatable

`T_a = 2π Aσ₀ / τ̇`. Locating the spectral peak requires `Aσ₀` in real stress
units. `ph-core::fault` currently emits stress **shape**, not magnitude — the
elastic response (Love numbers, radial structure) is not applied.

So even a perfect response spectrum could not be compared against the prediction
quantitatively. **F2 is no longer a parallel nicety; it is on the critical path.**

---

## Immediate — finish Phase 2 honestly

| # | Task | Why now |
|---|---|---|
| ~~C3b~~ | **Done.** 2,000 trials: slope 3.56, null median 0.43, null max 4.09, **p = 0.0095**. Survives, but weaker than the floor suggested | ✓ |
| **C2b** | Alias analysis — enumerate catalogue periodicities, compute beats against the constituent list, blacklist collisions | Outstanding validity gate. The 24 h detection cycle is strong enough (S1 power 16,245) that its beats will contaminate elsewhere. |
| ~~C4~~ | **Done** via analytic Doodson phases and a per-block shift null. 4/9 constituents survive BH-FDR. After normalising by forcing amplitude, **R(ω) is flat to within ~3× across 0.5–27 d** — the apparent band limitation was the tidal potential's own amplitude spectrum. No `T_a` located | ✓ |

C4 is the scientific centre of gravity here. It converts "LFEs respond to tides"
into "LFEs respond *like this* as a function of frequency", which is the object the
whole project is built to measure.

## Next — make the results transferable

| # | Task | Why |
|---|---|---|
| ~~F2~~ | **Done.** `ph_core::love`, degree-2 elastic calibration, good to ~2×. M2 solid Earth tide = **595 Pa** against Thomas et al.'s `Aσ₀` = 600 Pa, matching to 1% — independently explaining why Parkfield is non-linear and why C3's slope exceeds 2 | ✓ |
| ~~C5~~ | **Done.** 678,084 Cascadia tremor detections, 2009–2024. **8/9 constituents give the same verdict as Parkfield; M2, N2, O1 significant at both.** Different tectonics, geography and detection method, so a shared artifact is very hard to sustain. First replicated result in the project | ✓ |

## Then — Phase 3, which needs the split

| # | Task |
|---|---|
| **B1** | `pyo3`/`maturin` bindings on `ph-core` |
| **B2** | CLI emitting features as CSV/Parquet with provenance |
| **B3** | Create `EarthquakeForecastModule`, add `HANDOFF.md`, submodule `ph-core` |
| **P3.1** | USGS ComCat ingestion + magnitude-of-completeness per region and epoch |
| **P3.2** | GCMT focal mechanisms → real ΔCFS per event. **Now the highest-value step.** Raw tidal phase is blind to whether the tide loads or unloads each fault, so compressional and extensional responses cancel in the pooled statistic. Raises signal by aligning the feature with the physics rather than by discarding data |
| **P3.3** | ETAS baseline, fitted and frozen |
| **P3.4** | **Response spectrum for ordinary crust** — the band prediction test |
| **P3.5** | β(x,t) sensitivity field; CSEP evaluation |

**P3.4 is the whole point.** Everything before it is instrument.

---

## Critical path

```
C3b ─┐
C2b ─┼─> C4 ──┬─> C5 (independence)
     │        └─> F2 ──> P3.4  (band prediction, quantitative)
     │
     └────────> B1,B2,B3 ──> P3.1,P3.2,P3.3 ──> P3.4
```

**F2 and C4 are the two things that gate the central question.** Neither is large.

## Decision points

1. **After C4** — if the Parkfield response spectrum is *flat* across usable bands,
   the transfer-function framing is wrong and doc 08 §12 needs rethinking before
   Phase 3.
2. **After C5** — if Cascadia disagrees with Parkfield, the results are
   site-specific and the generalisation to earthquakes is unsupported.
3. **After F2** — if `Aσ₀` at Parkfield does not come out near Thomas et al.'s
   6×10⁻⁴ MPa, the stress pipeline has a scale error and every magnitude-dependent
   claim needs revisiting.
4. **Before P3.4** — pre-commit to what publishes if the band prediction fails. The
   honest answer is that a measured transfer function for ordinary crust is
   publishable *whatever shape it has*, and saying so now removes the incentive to
   keep slicing until something crosses p < 0.05.

## Standing rule, earned the hard way

**Specify the null before running, and state what structure it preserves.**

All four traps came from a null chosen after seeing the data. C2 was the first test
with the null fixed in advance, and it is the first result that survived
correction. That is not a coincidence.
