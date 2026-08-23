# 19 — Results

Consolidated findings. Written at the stopping point defined in
[16-plan.md](16-plan.md) decision 4, which pre-committed — before any result was
known — that the outcome would be reported whatever its shape.

Three separable contributions: a **methodological** one, a **positive** measurement,
and a **bounded null**.

---

## Summary

We built an instrument for measuring tidal stress from ephemerides and applied it
to four seismic catalogues spanning the Moon and Earth.

**Tidal Coulomb stress measurably modulates slow seismicity**, replicated at two
independent sites, surviving a test designed to destroy it. **We cannot detect it
in ordinary earthquakes**, bounded below 3.88% by five independent routes. And the
central theoretical prediction we set out to test turns out to be **unconstrained
by one to two orders of magnitude** by any catalogue we could obtain.

The most transferable result may be neither: **six documented ways to obtain a
confident wrong answer**, the worst reporting p = 10⁻⁸⁹ where the correct value
was 0.70.

---

## 1. Instrument and validation

`ph-core` (61 tests) computes tidal tensors from DE440 ephemerides, resolves
Coulomb stress on fault geometry, converts to Pa via degree-2 Love numbers, and
supplies analytic Doodson constituent phases. `eqf-analysis` (27 tests) adds
catalogue ingestion and the measurements.

**Validation ladder**, ordered by signal strength:

| Rung | Catalogue | Events | Role | Outcome |
|---|---|---|---|---|
| 1 | Apollo deep moonquakes | 6,954 | known answer | **5/5 periodicities to <0.21%** |
| 2 | Parkfield LFEs | 1,528,117 | strong effect | M2/N2/O1 significant |
| 2 | Cascadia tremor | 678,084 | independent site | M2/N2/O1 significant |
| 3 | ComCat M5.5+ / GCMT | 25,962 / 18,310 | the question | null, bounded |

**Elastic calibration, independently confirmed.** The M2 solid Earth tide computes
to **595 Pa**; Thomas et al. (2012) infer `Aσ₀ = 600 Pa` at Parkfield from tremor
triggering. Agreement to 1% from unrelated routes — and it explains why Parkfield
sits in the non-linear response regime (`S_T/Aσ₀ ≈ 1`).

**Null calibration.** The block-shift null returns uniform p-values (Kolmogorov-
Smirnov) on synthetic catalogues with no tidal signal, including Hawkes clustering
and a strong diurnal detection modulation. Power: ~92% at ε = 5% for 8,000 events,
scaling as √N. The power curve and the reported bounds were derived by different
routes and agree.

---

## 2. Positive result — tidal modulation of slow seismicity

### Replicated across two sites

| Constituent | Period (d) | Parkfield | Cascadia |
|---|---|---|---|
| M2 | 0.5175 | 223× null | 100× null |
| N2 | 0.5274 | 45.7× | 7.8× |
| O1 | 1.0758 | 137× | 51× |
| Mf–Sa | 13.7–365 | not significant | not significant |

**8/9 constituents give the same verdict.** The sites differ in tectonic setting
(strike-slip transform vs subduction megathrust), geography (~1000 km), epoch, and
**detection method** (template matching vs envelope cross-correlation).

### Amplitude scaling, faster than linear

Binning 1.53M events by ΔCFS cycle amplitude: response rises **71× across a 3×
amplitude range**, log-log slope **3.56** against 2 for linear. The trend survives
a re-binned null at **p = 0.0095**. Consistent with `R = R₀exp(S_T/Aσ₀)` in its
non-linear regime — which the independent elastic calibration predicts.

### The artifact hypothesis, tested and refuted

Ocean loading modulates microseism noise; Custodio et al. (2003) measured that
modulation **peaking at M2**. Detection threshold tracks noise, so detection
capability oscillates at our exact frequency — an explanation requiring zero
triggering, and fitting our pattern (signal in threshold-limited catalogues, none
in the complete one).

Stratifying by detection strength refuted it:

- **Parkfield by `ccsum`**: 22.81% (weakest) → 19.96% (strongest). **Ratio 0.87.**
  A detection artifact requires collapse toward zero.
- **Cascadia by magnitude**: 3.74% → 37.29%. **Ratio 9.98 — modulation rises
  tenfold with event size**, the opposite of the artifact prediction.

The magnitude dependence independently reproduces **Ide, Yabe & Tanaka (2016)**,
who report tides modulating the earthquake size distribution. Recovered from a
different catalogue with different machinery, while attempting to disprove
ourselves.

### Transfer function

Response per unit stress, with total tide (solid + ocean loading):

| Site | R(M2), per Pa |
|---|---|
| Parkfield | 4.45 × 10⁻⁴ |
| Cascadia | 2.33 – 3.05 × 10⁻⁴ |

**Agreement within a factor of two** across sites differing in geometry, setting,
detection and forcing amplitude — and robust to a 31% spread from assumed fault
geometry. `R` behaves like a property of the physics rather than of the site.

---

## 3. Bounded null — ordinary crust

Five independent routes, all null:

| Route | Result |
|---|---|
| Raw tidal phase, global, longitude-corrected | null |
| Depth stratification (pre-registered, Métivier prediction) | null; **prediction unsupported**, 2/6 bands |
| Mechanism-resolved ΔCFS sign (GCMT, both nodal planes) | null; 49.3% vs null median 49.2% |
| Total-tide R(ω) | null |
| Long-period bounds | unconstrained |

**Ordinary crust responds below 3.88% at M2 and 4.33% at O1**, where tremor shows
21.7% and 14.5% — **at least 3–5× weaker**.

Consistent with Beeler & Lockner's nucleation argument, and with why the literature
is genuinely mixed.

---

## 4. The central prediction is untestable with available data

Two timescales bound the responsive band: nucleation duration `t_n` (Beeler &
Lockner extrapolate ≥1 yr for the San Andreas) and Ader's critical period
`T_a = 2πAσ₀/τ̇` (~20–200 yr for ordinary crust). The prediction is that ordinary
crust responds at **years to decades** and is damped below.

Measuring R(ω) from M2 (0.52 d) to Ssa (183 d) with total tide:

| Site | R(M2) | long-period bound | ratio |
|---|---|---|---|
| Parkfield | 4.45e-4 | 3.17e-3 | 7× |
| Cascadia | 2.33e-4 | 4.24e-2 | 182× |

**The bound exceeds R(M2) by one to two orders of magnitude.** Neither refuted nor
confirmed.

The reason is structural: long-period forcing is **8–72 Pa** against **400–680 Pa**
at M2, and no long-period constituent is significant. Closing it requires
substantially more events (bounds scale as 1/√N) or stronger long-period forcing.
Neither exists in these catalogues.

---

## 5. Methodological findings

**Six ways to obtain a confident wrong answer**, each caught only because the result
was checked against what the statistic was actually testing:

| # | Failure | Naive answer | Correct |
|---|---|---|---|
| 1 | Time-shift null degenerate when catalogue shares the forcing's period | p = 10⁻⁸⁹ | **p = 0.70** |
| 2 | Uniform-time null tests temporal clustering, not alignment | 73/74 nests | **0/74 after FDR** |
| 3 | Raw period folding measures the detector | strong peaks | S1 artifact power 16,245 |
| 4 | Per-bin nulls when the binning variable derives from the forcing | highest bin p = 0.45 | trend nulled instead |
| 5 | Time-shift null against a demodulated constituent — **trap 1 recurring, walked into after documenting it** | — | block-shift required |
| 6 | Sham-frequency null returns the leaking constituent's phase | D²/N ≈ N where no tide exists | not a baseline |

**All six share one shape: a statistic that looks decisive while silently answering
a different question.**

Two further failures were not nulls at all:

- **A raw response measurement is not a transfer function.** Dividing by the forcing
  is the definition, not a refinement. Skipping it produced a confident,
  physically-flavoured, wrong conclusion that survived a full write-up.
- **Agreement between methods is not agreement about physics when the methods share
  an input.** Three independent nulls agreed while sharing the same incomplete
  forcing.

**This may explain the literature's mixed record.** If null choice can move an
answer by 88 orders of magnitude, studies differing only in null construction will
disagree — and will each look internally sound.

---

## 6. What is unresolved

| Item | Status |
|---|---|
| Band prediction | **Unconstrained**; needs more events or stronger long-period forcing |
| Cross-band loading phase calibration | Verified for **M2 only**; `hartid`'s constituent inference blocks validation |
| Cascadia fault geometry | **Assumed.** Not load-bearing for R(ω), but is for the lag comparison |
| Semidiurnal/diurnal lag split (132.6°) | Real, coherent (diurnal agrees to 1.6°), **unexplained**. Ocean loading refuted as the cause |
| External validation | Thin. Ide reproduction is one genuine instance; Métivier remains untried |

---

## 7. Reproducibility

All data is public and requires no credentials. Every dataset has a fetch script.

```bash
git submodule update --init --recursive
./scripts/fetch-kernels.sh && ./scripts/fetch-apollo.sh
./scripts/fetch-parkfield.sh && ./scripts/fetch-cascadia.sh
./scripts/fetch-comcat.sh && ./scripts/fetch-gcmt.sh
./scripts/setup-spotl.sh          # ocean loading
cargo test && cargo run --release --example moonquake_periodogram
```

Sources: NASA PDS (Apollo), USGS ScienceBase (Parkfield LFEs), PNSN (Cascadia
tremor), USGS ComCat, Global CMT, NAIF (DE440), SPOTL with GOT4.7 and FES2004.

Two API behaviours are documented in the fetch scripts because each would silently
corrupt a catalogue: the PNSN service **caps responses at 20,000 events with HTTP
200 and no truncation flag**, and returns 404 rather than an empty result for
windows with no data.


---

## Sequel: the astrological programme

The work above tests physically motivated tidal hypotheses. A separate programme
then tested ~9,800 astrological chart features against the same catalogues, with
no physical filtering and no requirement that any mechanism be understood.

Results in [23-ml-results.md](23-ml-results.md) and
[25-stratified-results.md](25-stratified-results.md); pre-registered in
[21](21-preregistration.md) and [24](24-stratified-preregistration.md).

It found nothing, and the bounds are quantitative: **~4% per standard deviation
globally, ~6–9% within any focal-mechanism class.** The global bound is
comparable to the 3.88% M2 bound reached above by five independent physical
routes, which is a useful check — two entirely different approaches converging on
the same order of magnitude.

Two findings from that programme bear on the physics work here:

- **The effective sample size of global seismicity is ~10,000 independent events
  over 49 years, not 488,000.** Windowed declustering removes aftershock
  sequences but leaves long-range regional correlation; forcing 500 km and 365
  days of separation is what it takes to calibrate a permutation null. Any study
  quoting statistics on hundreds of thousands of events is quoting a sample size
  it does not have.
- **Annual counts are over-dispersed 2.4× even at M6.0+**, so any analysis
  comparing across years — including much of the long-period tidal literature —
  is confounded with the rate history unless it models it.
