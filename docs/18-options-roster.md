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

---

# Expansion — 2026-08-22

More options, and more detail on how the leading ones would actually be done.

---

## §G — Natural experiments

A coherent theme that has never been exploited. Rate-and-state gives
`T_a = 2π Aσ₀ / τ̇`, so **when stressing rate rises, `T_a` falls and tidal
sensitivity should rise.** Several settings vary `τ̇` by orders of magnitude *within
our existing data*.

| # | Item | Effort | The prediction |
|---|---|---|---|
| **G1** | **Aftershocks vs background** | M | During an aftershock sequence the fault sits at criticality with elevated `τ̇`. Tidal sensitivity should be **higher in aftershocks than in background events.** Declustering already needed for ETAS; this reuses it. Sharp, cheap, and uses data in hand. |
| **G2** | **Inside vs outside ETS episodes** | M | Cascadia tremor occurs in episodic tremor-and-slip bursts every ~14 months, during which local stressing rate is enormously elevated. Cleaner separation than G1 — episodes are unambiguous in the catalogue. |
| **G3** | Injection-induced seismicity | L | Oklahoma, Groningen: **known, recorded stressing history.** Identified in doc 05 §3 and never used. The closest thing to a controlled experiment available. Also a methodological check — if the pipeline attributes injection-driven seismicity to tides, that is a decisive failure caught cheaply. |
| **G4** | Volcanic / geothermal seismicity | M | High pore pressure, compliant media, expected high sensitivity. Another point on the `Aσ₀` axis. |

**G1 and G2 are the strongest untested predictions in the project** after the band
prediction, and unlike it they are answerable with data already downloaded.

## §H — Physical consistency checks

Cheap tests that a real effect must pass and an artifact need not. None has been run.

| # | Item | Effort | Why |
|---|---|---|---|
| **H1** | **Cross-constituent phase consistency** | S | If triggering is real, preferred phase at M2, N2 and O1 should reflect a *single* physical lag, not three unrelated numbers. Scattered phases would indicate something is wrong. We have all three phases already. |
| **H2** | **Stress *rate* vs stress amplitude** | S | Rate-and-state says `dΔCFS/dt` matters alongside `ΔCFS`; Weber found rates mattered for some moonquake clusters. We have only ever tested amplitude. A distinct pre-registered alternative, and the CLI already emits `dcfs_dt`. |
| **H3** | Amplitude law across sites | S | Parkfield shows ε = 21.7%, Cascadia 12.7%. If `R = ε/stress` is a property of the physics, the difference should be explained by the stress amplitudes at each site. A two-point transfer-function check across sites. |
| **H4** | Perigee–apogee amplitude knob | S | Lunar distance varies 5.5%; tides go as 1/d³, giving **18% amplitude modulation** at 27.55 d with known phase. Designed in doc 08 §13c, never run. |
| **H5** | 18.61 yr nodal envelope | M | Modulates diurnal amplitudes ±11%, giving a slow envelope no instrumental effect plausibly mimics. Parkfield's 23 years gives only ~1.2 cycles — marginal. Better on a longer catalogue. |

**H1 and H2 are hours of work on data already in memory.**

## §I — A novel application, independent of the triggering question

| # | Item | Effort | Why |
|---|---|---|---|
| **I1** | **Tidal ΔCFS as a nodal-plane discriminator** | M | A moment tensor gives two planes and does not say which broke. If tidal triggering is real, **the plane showing stronger ΔCFS preference is more likely the true fault.** Weber et al. did exactly this for moonquakes. For earthquakes it can be *validated*: plenty of events have the plane independently determined from aftershock distributions or surface rupture. |

This is worth flagging separately because it **inverts the problem**. Instead of
using known geometry to test triggering, it uses triggering to infer geometry — and
it is falsifiable against ground truth, which the triggering question itself is not.
It would also be useful to seismology regardless of how the band prediction turns
out.

---

## How the leading items would actually be done

### A1 — R(ω) with ocean loading

1. Run `ocean-loading-sites.sh` at **Parkfield and Cascadia** for M2, N2, O1, Q1,
   Mf, Msf, Mm, Ssa, Sa. Two sites, nine constituents — minutes, not the 18,316-site
   job.
2. Reconstruct loading strain at event times:
   `ε(t) = A cos(χ_local(t) + φ_SPOTL + 10.97°)`, verified in
   `verify_loading_phase.rs`.
3. Strain → stress under the free surface (`σ_zz ≈ 0`, plane stress with
   `μ = 30 GPa`, `ν = 0.25`), giving the horizontal tensor.
4. Add to the solid-tide tensor, resolve on fault geometry.
5. Per constituent, fit total forcing amplitude by least squares on the analytic
   argument — the same method already used for solid tide.
6. `R(ω) = ε(ω) / amplitude_total(ω)`. Compare against the solid-only version.

**The decisive comparison:** if R(ω) was flat with solid tide and *rises* with total
tide, the band prediction was killed by an artifact of our own making.

### B1 — Reproducing Métivier et al. (2009)

Their result: NEIC, 442,412 events, ~99% confidence, events preferentially at
**ground uplift** — reduced normal stress. Anomaly larger for smaller and shallower
events.

1. NEIC is ComCat's source, so the catalogue is already reachable. Match their
   magnitude and epoch cuts.
2. Compute solid-tide **normal stress** at each event (not Coulomb — they used the
   tidal potential's vertical component).
3. Schuster on the phase, plus our block-shift null for a modern comparison.
4. Check three things: the **sign** (uplift), the **confidence**, and the reported
   **depth and magnitude dependence**.

**Why it matters more than another internal test:** it validates the *stress* path
end to end against someone else's published number. The moonquake test validated
timing only. Their reported depth dependence is also a second, independent check on
our P3.5 null, which found no depth dependence.

### A4 — The 10.97° diagnostic

10.97° of M2 is 22.7 minutes. If it is a **constant time offset**, O1 (25.82 h)
must show `22.7/1549 × 360 = 5.3°`. If it is a **phase convention**, O1 shows ~11°.

The O1 loading data is already computed. One least-squares fit against `hartid`
output separates them. A time offset would mean a real timing bug somewhere; a
phase convention is benign.

---

## §J — Library, expanded

| # | Item | Effort | Why |
|---|---|---|---|
| **J1** | Doodson generalised to arbitrary body pairs | M | The original "not just base 12" idea. Extend from the six tidal arguments to arbitrary integer combinations of planetary longitudes under the d'Alembert constraint (doc 02). `doodson` already has the machinery; this is the astrology-facing generalisation. |
| **J2** | Golden reference dataset + regression tests | S | The library has no protection against silent numerical drift. Freeze known-good outputs. |
| **J3** | Columnar batch benchmarks across the WASM boundary | S | Doc 06 asserts boundary crossings dominate. Never measured. |
| **J4** | Determinism harness | S | Seeds, kernel versions, commit hashes recorded with every result, so any number can be regenerated exactly. |

---

## Revised recommendation

The roster is now large enough that ordering matters more than completeness.

```
A4          minutes, diagnostic, unblocks confidence in timing
 ├─> A2 -> A1 -> A3        the live scientific question
 ├─> H1, H2                hours; consistency checks that should already exist
 ├─> G1, G2                strongest untested predictions, data in hand
 └─> B1                    external validation, independent of everything else
```

**A4, H1, H2 in one session** — all small, all diagnostic, and H1/H2 could
independently undermine or strengthen the tremor result before more is built on it.

**Then A1**, because every interpretation currently rests on a control we know is
compromised.

**Then G1/G2**, because they are the strongest untested predictions we can actually
answer, and unlike the band prediction they do not need a larger catalogue.

**B1 in parallel**, because it is the only item that checks this work against
someone else's.
