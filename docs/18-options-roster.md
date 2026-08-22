# 18 — Options Roster

Written 2026-08-22, after the self-critique closed two of its own points. This is
the full option space, not a plan — pick from it.

**Effort:** S = hours, M = a session, L = several sessions.

---

## The fork worth deciding first

Two coherent programmes, and drifting between them serves neither.

| | **Seismology** | **Library** |
|---|---|---|
| Question | Does tidal stress modulate seismicity? | Bulk harmonic computation for four apps |
| Status | One replicated positive, one bounded null, one open question | Physics core built; the stated deliverable is not |
| Serves | EarthquakeForecastModule | Star Seer, Resonant Finder, Cosmic Cypher, AstrologyCore |
| Next | §A, §B | §D |

The original brief was the library. Nearly all effort has gone to seismology. That
is not necessarily wrong — the seismology produced `doodson`, which is the core
Star Seer primitive, arrived at by accident — but it should be a choice.

---

## §A — Finish what is in flight (seismology)

| # | Item | Effort | Why |
|---|---|---|---|
| **A1** | **Recompute R(ω) at tremor sites with ocean loading** | M | **The live question.** Loading is large at semidiurnal and negligible at long period, so correcting it lowers R at short periods only — turning flat into *rising*, which is what the band prediction predicts. The current flat-R(ω) claim is circular and must be retracted or confirmed. |
| **A2** | Ocean loading for Mf, Mm, Ssa, Sa | S | A1 is incomplete without the long-period end. Data generation is ~6 min per constituent. |
| **A3** | Total-tide earthquake test (P3.2 redo) | M | Determines whether three nulls survive contact with the omitted forcing. Meaningless before A1, since it needs a control we trust. |
| **A4** | Diagnose the 10.97° offset | **S** | Constant *phase* or constant *time*? O1 would show 5.3° if it is a time offset, ~11° if a phase convention. Data already computed. Tells us whether a timing bug lurks elsewhere. |
| **A5** | Bound Parkfield's residual detection effect | S | The 0.87 ccsum ratio is consistent with a small detection component riding on real triggering. Quantify rather than dismiss. |

**A4 first** — minutes, and it is diagnostic. Then A2 → A1 → A3.

## §B — New science

| # | Item | Effort | Why |
|---|---|---|---|
| **B1** | **Reproduce Métivier et al. (2009)** | M | **The external validation we lack.** NEIC, 442k events, published ~99% confidence and phase preference toward uplift. Validates the *stress* path, which the moonquake test never touched. Every check so far is my code against my code. |
| **B2** | β(x,t) time-varying sensitivity | L | The strongest untested idea in the project, and now testable — both tremor sites have strong signal, and Beaucé reports sensitivity rising 1.5 yr before Ridgecrest. This is the part of the field with forecasting value. Pre-register window length and excursion criterion. |
| **B3** | Magnitude-distribution (b-value) modulation | M | We found 10× size-dependence at Cascadia while refuting detection bias. Ide et al. predict tides modulate the size distribution. Do it deliberately rather than as a by-product. |
| **B4** | Third site: Nankai / Japan tremor | M | Parkfield and Cascadia are both west-coast North America, both analysed with my code. A Japanese catalogue is genuinely independent in network, operator and region. |
| **B5** | Spectral reformulation — Lomb-Scargle with red-noise null, or PTA Gaussian process | L | **Methodologically superior frame.** Recasts "do events cluster at phase φ(t)" as "does this point process have excess power at ω beyond its own red spectrum". Dissolves the whole trap class rather than patching it, and never assumes independence. I reached for invention when the field had the tool. |
| **B6** | Deep moonquake Coulomb against Weber's published planes | M | Phase 1 left 0/74 surviving FDR with a 7.1σ ensemble excess. Weber's per-cluster constraints are a second external check. |

## §C — Methodology debt

| # | Item | Effort | Why |
|---|---|---|---|
| **C1** | Project-wide multiple-comparison accounting | S | FDR was applied within tests, never across the programme. Many tests have been run. |
| **C2** | Fix FDR over non-independent tests | S | "9/12 families survive FDR" is invalid as applied — the families are co-located and not independent. Either use a dependence-aware procedure or restate the claim. |
| **C3** | Publish a power curve with every null | S | The calibration run showed bounds and power agree; make that routine rather than incidental. |
| **C4** | Pre-registration document for remaining tests | S | Two of the last three tests were pre-registered and both were informative. Formalise it. |

## §D — The library (upstream, PlanetaryHarmonicsModule)

The stated deliverable. None of this exists.

| # | Item | Effort | Why |
|---|---|---|---|
| **D1** | **WASM / TypeScript surface** | L | The original brief. Rust-to-Rust composition into one WASM module is already supported by RustSPICE; nothing blocks it. |
| **D2** | Angle-domain root finding | M | Star Seer's millisecond micro-aspect events. Constituent arguments are near-linear in time, so crossings solve analytically and Newton-refine — O(events) instead of O(sample rate). `doodson` already provides the arguments. |
| **D3** | Harmonic ephemeris precompute | M | O(1) timestream queries from a precomputed (frequency, phase, amplitude) table. Resonant Finder needs this. |
| **D4** | HEALPix global synthesis | M | Degree-2 fields are 5 coefficients; global evaluation should be O(coefficients), not O(locations). Deferred once already. |
| **D5** | Multi-basis harmonic encoding | M | The original "not just base 12" idea. Fourier features, d'Alembert constraint, group-lasso selection — designed in doc 02, never built. |
| **D6** | Move ocean loading into `ph-core` | M | It is physics and it is shared. Currently a shell script against vendored Fortran. |

**D2 is the sleeper.** The seismology work produced `doodson` — analytic angular
arguments, validated to 0.1% against published constituent periods, with the
longitude correction independently confirmed at four sites. That is exactly Star
Seer's core primitive, already built and tested.

## §E — Engineering

| # | Item | Effort |
|---|---|---|
| **E1** | Vendor SPOTL reproducibly (currently a scratchpad build) | S |
| **E2** | CI on both repos | S |
| **E3** | Data checksums and provenance manifests | S |
| **E4** | Benchmarks for the bulk paths | S |

## §F — Output

| # | Item | Effort | Why |
|---|---|---|---|
| **F1** | **Methodological write-up** | M | Six traps, each producing a confident wrong answer — the worst p = 10⁻⁸⁹ against a true 0.70 — plus a calibrated null and a refuted artifact hypothesis. **This may be the most transferable thing the project has produced**, and it explains why the tidal-triggering literature is mixed. |
| **F2** | Tremor result write-up | M | Two sites, magnitude dependence reproducing Ide et al., detection bias excluded by test. |
| **F3** | Upper-bound write-up | M | Pre-committed as publishable whatever the result. **Blocked on A1/A3** — the bound is currently solid-tide-only. |

---

## Recommended sequence

```
A4  (minutes, diagnostic)
 └─> A2 -> A1 -> A3        the live scientific question
      ‖
      B1                    external validation, independent of A
      ‖
      D2 -> D1              library, and D2 is nearly free given doodson
```

**If forced to pick three:** A4, A1, B1. The first is cheap and diagnostic, the
second is the only live question, the third addresses the deepest remaining
weakness — that nothing here has been checked against anyone else's result.

**If the goal is the original brief instead:** D2, D1, D5. The seismology has
already produced the hard part of D2.
