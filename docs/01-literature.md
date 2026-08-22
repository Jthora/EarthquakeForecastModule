# 01 — Tidal Triggering Literature

The honest summary: **tidal triggering of earthquakes is real, small, and strongly
dependent on tectonic setting.** Effect sizes range from undetectable (Southern
California, general catalogues) to a factor-of-3 rate change (shallow thrust
faults under high tidal stress). Anything claiming a large universal effect is
wrong; anything claiming zero effect is also wrong.

---

## Foundational positive results

### Tanaka, Ohtake & Sato (2002), *JGR* 107(B10)
*Evidence for tidal triggering of earthquakes as revealed from statistical analysis of global data.*

Global catalogue analysis using the **Schuster test** on tidal phase angles.
Established the standard methodology the field still uses. Found significant
phase selectivity, strongest in the period preceding large events.

**Why it matters here:** the Schuster test is the field's canonical statistic and
is exactly the *n = 1* case of the Fourier angular encoding in
[02-angular-encoding.md](02-angular-encoding.md). Our generalisation to arbitrary
harmonic order *n* is a direct, defensible extension reviewers already understand.

### Cochran, Vidale & Tanaka (2004), *Science* 306, 1164
*Earth tides can trigger shallow thrust fault earthquakes.*

The strongest positive result in the literature. Shallow thrust faults in
subduction zones, where ocean tidal loading produces large stress amplitudes.

- Earthquake rate varied **by a factor of ~3** with tidal stress
- Probability of chance occurrence **< 1 in 10,000**
- Best correlation at friction coefficient **μ = 0.4** (good for μ = 0.2–0.6)

**Key lesson:** the effect is large *where tidal stress amplitude is large*. Setting
selection is not cherry-picking here — it is the physical prediction. Ocean
tidal loading, not solid Earth tide, drives this.

### Métivier et al. (2009), *EPSL* 278, 370–375
*Evidence of earthquake triggering by the solid earth tides.*

Largest global catalogue analysis: **NEIC, 442,412 events**.

- Clear correlation at **~99% confidence**
- Events occur slightly more often during **ground uplift** (reduced normal stress)
- Anomaly is **larger for smaller and shallower** earthquakes

**Key lesson:** confirms the sign of the mechanism (normal stress unclamping) and
shows the effect concentrates in small, shallow events — where catalogues are
largest but completeness is worst. Magnitude-of-completeness control is critical.

### Scholz, Tan & Albino (2019), *Nature Communications* 10, 2526
*The mechanism of tidal triggering of earthquakes at mid-ocean ridges.*

Mid-ocean ridges show strong triggering, but with a **surprising phase**: events
peak at low tide, opposite to the naive expectation. Explained by magma chamber
inflation/deflation rather than direct fault stress.

**Key lesson:** phase, not just amplitude, encodes mechanism — and the correct
phase can be counter-intuitive. This is a direct argument for keeping the sine
component in angular encodings so a model can *learn* phase rather than assume it.

---

## The mechanism constraint

### Beeler & Lockner (2003), *JGR* 108(B8)
*Why earthquakes correlate weakly with the solid Earth tides: effects of periodic stress on the rate and probability of earthquake occurrence.*

The most important methodological paper for this project. Laboratory faults
loaded by constant stressing rate plus a small sinusoid, interpreted through
rate-and-state friction.

Findings:
- Correlation is facilitated only when the **forcing period exceeds the
  characteristic nucleation time** of frictional instability. Tidal periods are
  near or below it, hence the weak correlation.
- **~10⁴ earthquakes are required** to demonstrate a statistically robust
  correlation. ⚠ *Corrected in the fourth pass — this doc originally said 10⁵–10⁶,
  taken from a search summary and never verified.* The paper's equation 18 gives
  `N ≥ ln(P_rw)/(Δτ_u/(2aσ_n))²`; worked examples give **6.2×10³–5.5×10⁴** events
  for Δτ = 0.01 MPa, and the abstract states daily Earth tides require
  **">13,000 earthquakes to detect."**
- A **10× increase in tidal stress amplitude gives a 100× decrease** in the number
  of events needed for detection.

**Direct consequences for this project:**

1. The full USGS catalogue pull remains worthwhile — but for **coverage and
   confounder control, not raw statistical power**. At ~10⁴ events needed, ComCat's
   millions clear the bar easily and many *regional* subsets are viable too. Re-run
   the doc 04 §7 power analysis with the inverse-square law.
2. Prioritise **high-amplitude-stress settings** (subduction, coastal, ridges)
   where the quadratic amplitude advantage collapses the sample requirement.
3. The nucleation-time argument predicts **long-period constituents should
   correlate better than semidiurnal ones** — and the fourth pass sharpened this
   dramatically. Beeler & Lockner give the nucleation duration **t_n ≥ 1 year for
   the San Andreas**, damping response above 1/t_n. Combined with Ader's critical
   period T_a ~ 20–200 yr, the predicted responsive band for ordinary earthquakes
   is **~1 year to ~200 years** — Sa at the edge, the 18.61 yr nodal term and
   decadal LOD inside, and **all semidiurnal, diurnal, fortnightly and monthly
   constituents damped.** See doc 07, fourth pass.

Point 3 is the single most useful physical hook in this literature. It gives a
mechanism-grounded reason to look at long-period orbital terms.

---

## Stress-state proxy — the strongest reframe

### Beaucé et al. (2023), *GRL* 50, e2023GL104375
*Enhanced tidal sensitivity of seismicity before the 2019 M7.1 Ridgecrest earthquake.*

Built a 10-year, **>150,000 event** catalogue preceding Ridgecrest. Found a robust
increase in tidal sensitivity of seismicity along the fault beginning
**~1.5 years before** the mainshock.

### Ide, Yabe & Tanaka (2016), *Nature Geoscience* 9, 834–837
*Earthquake potential revealed by tidal influence on earthquake size–frequency statistics.*

Very large earthquakes (2004 Sumatra, 2010 Maule, 2011 Tohoku-Oki) tend to occur
near times of **maximum tidal stress amplitude**. The tendency is not obvious for
small events — i.e. tides modulate the **size distribution (b-value)**, not just
the rate.

**Why these two reframe the whole project:**

The naive target is "do tides trigger earthquakes?" — a small effect requiring
enormous samples. The better target is:

> **Tidal sensitivity is a measurable proxy for how close a fault is to failure.**

A fault near critical stress responds more to a given small perturbation. So the
*time-varying responsiveness* to tidal forcing carries information about stress
state — and that is a precursor signal, not merely a triggering correlation.

This changes the ML target from "rate modulation" to "**time-varying coupling
coefficient between tidal stress and seismicity rate**," which is a far richer and
more forecast-relevant quantity. See [04-ml-architecture.md](04-ml-architecture.md).

---

## Negative and contested results

Include these. A literature base that only cites positive results is not a
literature base.

- **Southern California** (*GJI* 205, 681): no significant tide–earthquake
  correlation for any studied catalogue. A well-studied, well-instrumented region
  with a null result.
- **Low-frequency earthquakes near Parkfield** (Thomas et al. 2012, *JGR* 117):
  some apparently significant tidal correlations may be **spurious**, driven by
  correlation with other stress components rather than tidal triggering.
- **Bendick & Bilham (2017)**, *GRL* 44, 8320–8327, *Do weak global stresses
  synchronize earthquakes?* — reports a ~32-year cycle in M7+ rate correlated with
  **length-of-day** variations, with a ~5-year lead. Published and interesting;
  widely questioned on the small number of effective cycles in a 115-year record.
  Treat as Tier C.
- **Holt & Newman (2025)**, arXiv:2508.07064 — Jupiter/M7+ claim. See the
  multiple-comparisons analysis in [00-framing.md](00-framing.md). Cite only as a
  methodological counter-example.

---

## Reference table

| Work | Setting | Effect | Tier |
|---|---|---|---|
| Tanaka et al. 2002 | Global | Phase selectivity, Schuster test | A |
| Cochran et al. 2004 | Shallow thrust / subduction | **Rate ×3**, p < 10⁻⁴ | A |
| Métivier et al. 2009 | Global, 442k events | ~99% conf., uplift phase | A |
| Beeler & Lockner 2003 | Lab / rate-and-state | Needs 10⁵–10⁶ events | A |
| Scholz et al. 2019 | Mid-ocean ridges | Strong, inverted phase | A |
| Ide et al. 2016 | Global, M8+ | b-value modulation | A/B |
| Beaucé et al. 2023 | Ridgecrest | Sensitivity rises 1.5 yr prior | A/B |
| Rubinstein et al. 2008 | Cascadia tremor | **Very strong** (see doc 05) | A |
| S. California (GJI 205) | Regional | Null | A |
| Bendick & Bilham 2017 | Global, LOD | 32-yr cycle, contested | C |
| Holt & Newman 2025 | Global, Jupiter | Multiple-comparison artifact | — |

## Open questions worth owning

1. Does the Beeler–Lockner nucleation argument actually predict better correlation
   at long periods (Mf, Mm, Ssa, Sa)? Testable directly with our feature basis.
2. Can tidal sensitivity be estimated continuously in space and time, rather than
   per-region retrospectively? This is the Beaucé result generalised, and it is
   an ML problem.
3. Is the b-value modulation of Ide et al. recoverable as a learned magnitude
   distribution conditioned on tidal amplitude?
