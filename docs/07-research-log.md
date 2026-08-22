# 07 — Research Log

Dated record of literature passes, findings, and decisions. Append-only.

---

## 2026-08-21 — First literature pass

**Scope:** tidal triggering of earthquakes; angular feature encoding; ML
forecasting architecture; validation datasets.

### Headline findings

**1. The effect is real, small, and setting-dependent.**
Range spans null results (Southern California, general catalogues) to a
**factor-of-3 rate change** (Cochran et al. 2004, shallow thrust faults under high
ocean-tidal loading, p < 10⁻⁴). Métivier et al. 2009 found ~99% confidence on
442,412 NEIC events, with the anomaly *larger for smaller and shallower* events.

Neither "tides cause earthquakes" nor "tides have no effect" is supportable. The
defensible position is that the effect exists, is a few percent globally, and
grows large where tidal stress amplitude is large.

**2. Sample-size requirement is a hard design constraint.**
Beeler & Lockner (2003): **10⁵–10⁶ earthquakes** needed for a statistically robust
correlation. A 10× stress amplitude increase cuts the requirement 100×.

*Decision:* full USGS catalogue pull is a hard requirement. Power analysis before
any modelling. Region selection driven by tidal stress amplitude, which is a
physical prediction rather than cherry-picking.

**3. The nucleation-time argument predicts long periods should win.**
Beeler & Lockner: correlation is facilitated when forcing period **exceeds the
characteristic nucleation time** of frictional instability. Tidal semidiurnal
periods sit near or below it — hence weak correlation.

This predicts that **long-period constituents (Mf 13.66 d, Mm 27.55 d, Ssa 182.6 d,
Sa 365.26 d, nodal 18.61 yr) should correlate better than semidiurnal ones.**
Sharp, testable, and it is exactly the regime where orbital/planetary cycles live.

*This is the most useful physical hook found in the entire pass.* See
[08-hypotheses.md](08-hypotheses.md) §1 for what it becomes.

**4. Phase is diagnostic, and can be counter-intuitive.**
Métivier: events cluster at the **uplift** phase. Scholz et al. 2019: mid-ocean
ridge events peak at **low tide** — opposite the naive expectation, explained via
magma chamber inflation.

*Decision:* angular encodings must retain the sine component so phase is learned,
never assumed. Settles the encoding question in [02-angular-encoding.md](02-angular-encoding.md).

**5. Tidal sensitivity is a stress-state proxy — the strongest reframe available.**
Beaucé et al. 2023: tidal sensitivity of seismicity rose **~1.5 years before** the
M7.1 Ridgecrest mainshock (>150,000 event catalogue). Ide et al. 2016: very large
events cluster near maximum tidal amplitude, modulating the **size distribution**
rather than the rate.

*Decision:* primary ML target changed from rate modulation to the **time-varying
sensitivity field β(x,t)**. See [04-ml-architecture.md](04-ml-architecture.md) §4.

**6. Tremor is an order of magnitude easier than earthquakes.**
Rubinstein et al. 2008: clear tremor pulsing at 12.4 h and 24–25 h; tides
"influence the genesis of tremor much more effectively than they do the genesis of
normal earthquakes." Ide et al.: tremor rate rises **exponentially** with tidal
stress.

**7. Deep moonquakes are tidally dominated, not tidally nudged.**
Apollo PSE catalogue. Known monthly / 206-day / 6-year periodicities.

*Decision:* adopt a graded validation ladder — **moonquakes → tremor →
earthquakes**. Each phase has a known-answer test before the next. See
[05-research-frontier.md](05-research-frontier.md).

### Methodological finding

**Holt & Newman (2025), arXiv:2508.07064**, *Tidal Triggering of M7+ Earthquakes by
Jupiter*: **28,782 chi-squared tests**, **1,071 significant** reported. At α = 0.05
chance predicts ≈ 1,439. They found *fewer significant results than the null
expects*.

(Tests are heavily correlated — overlapping 5–45 day intervals — so the effective
independent count is well below 28,782 and 1,439 is not the exact expectation. The
conclusion is unchanged: uncorrected and unadjusted for correlation, the result is
indistinguishable from noise.)

*Decision:* FDR correction mandatory on any feature scan; time-shifted null
distribution rather than random shuffling. Documented as the project's canonical
cautionary example in [00-framing.md](00-framing.md).

### Physics check performed

Tidal effect scales as M/d³. Computed:

| Body | M/d³ (SI) | vs. Moon |
|---|---|---|
| Moon | 1.30 × 10⁻³ | 1 |
| Sun | 5.9 × 10⁻⁴ | 0.46 |
| Venus (closest) | 8.9 × 10⁻⁸ | 7 × 10⁻⁵ |
| Jupiter (closest) | 9.4 × 10⁻⁹ | 7 × 10⁻⁶ |

Matches published figures (Venus ~10⁴× weaker than Sun+Moon combined).

*Conclusion:* direct planetary tidal triggering is below detection and must be
documented as such. Planetary influence enters legitimately via **orbital
perturbation of the lunisolar tide** and possibly via **principal-axis deflection**
(see [03-tidal-tensor.md](03-tidal-tensor.md) §5).

*Also recorded:* a uniform gravitational field produces no stress. Barycentres,
Lagrange points, and inter-planetary "gravitational nulls" are not stress
mechanisms. Written into doc 03 explicitly to prevent the reasoning error.

### Tooling identified

| Tool | Purpose |
|---|---|
| `pyCSEP` (cseptesting.org) | N/S/M/L tests, information gain per earthquake |
| RECAST (`github.com/keliankaz/recast`) | Neural temporal point process, GRU encoder-decoder |
| EarthquakeNPP (arXiv:2410.08226) | Benchmark suite for neural point processes |
| SPOTL (Agnew) | Ocean tidal loading |
| GCMT | Focal mechanisms → fault geometry |
| IERS | Polar motion, LOD series |
| Apollo PSE (Nakamura) | Deep moonquake catalogue |

### Confounder identified — high priority

**Seasonal hydrological loading is annual. So are Sa and Ssa.** Any annual-period
"celestial" correlation is confounded with groundwater, snow load, and atmospheric
pressure by default. Must enter as an explicit covariate, not be ignored.

Currently the most dangerous unaddressed confounder in the project.

### Open threads for next pass

1. Rate-and-state transfer function — analytic form of amplitude and phase
   response to periodic forcing (Ader & Avouac; Dieterich 1994; Beeler & Lockner)
2. Cartwright–Tayler–Edden tidal potential catalogue — amplitudes per constituent
3. Farrell (1972) load Love numbers and Green's functions for ocean loading
4. Schuster test null distribution generalised to harmonic order n
5. ETAS space–time–magnitude likelihood, exact form
6. Deep moonquake catalogue structure and prior periodicity analyses

Carried into [09-deep-dive-agenda.md](09-deep-dive-agenda.md).

---

## 2026-08-21 — Second pass: rate-and-state response theory

**Scope:** agenda items 1 (transfer function), 2 (spectral prediction test),
PTA statistical methods. Triggered by hypothesis 9.

### Finding 1 — the spectral prediction in hypothesis 9 was wrong ✗

The first-pass guess was a simple low-pass with response ∝ 1/ω across the tidal
band, giving a monotonic ranking of constituents by period (nodal > Sa > Ssa > Mm >
Mf > semidiurnal).

**The literature does not support this.** The actual result (Ader, Lapusta, Avouac
& Ampuero 2014, *GJI* 198, 385–413; Heimisson & Avouac 2020, *GRL*) is a
**two-regime response** separated by a critical period.

Define the Dieterich characteristic time and critical period:

```
t_a = a σ / τ̇          T_c ≈ 2π t_a
```

- **T ≪ t_a** — seismicity rate **tracks stress**. Response amplitude ∝ A,
  approximately flat in period.
- **T ≫ t_a** — seismicity rate **tracks stressing rate**. Response amplitude
  ∝ A/T, *falling* with increasing period.

So long-period constituents are **not** automatically favoured. The nodal-term
ranking prediction is dead.

**What survives, and is sharper:** the response has a **knee, and plausibly a peak,
near T ≈ t_a**. A 2012 AGU abstract from the same group describes an amplitude
maximum at a characteristic period *decreasing at both smaller and larger periods*
— i.e. **band-pass**.

*Reconstruction (unverified — see caveat):* band-pass behaviour requires **two**
timescales, not one. The Dieterich relaxation time t_a governs the long-period
roll-off; a separate **nucleation duration** governs the short-period cutoff
(Beeler & Lockner's original argument). Response peaks between them.

If correct, this is *better* for the inverse problem than the original guess: the
spectrum would identify **two** fault properties rather than one, and a peak is far
more diagnostic than a slope — an artifact can produce a slope, but is unlikely to
produce a peak at a physically predicted location.

**Caveat, important:** Heimisson & Avouac state their analytical model is
**not valid for periods similar to t_a** — precisely the region of interest. The
band-pass reconstruction above is inference from abstracts, not from the papers
themselves. Both primary sources are paywalled and were not obtained this pass.
**Do not build on this until the actual equations are in hand.**

### Finding 2 — hypothesis 1 is already an active research programme

Frequency-dependent tidal sensitivity as a probe of fault properties is
**published work**, not a novel idea:

- **Beeler et al. (2018)**, *JGR* 123 — *Constraints on friction, dilatancy,
  diffusivity, and effective stress from low-frequency earthquake rates on the deep
  San Andreas Fault*
- **EPSL (2025)**, *Probing lower-crustal fault properties with frequency-dependent
  tidal tremor triggering* — San Andreas LFEs. Finds that **modulation amplitude is
  controlled mainly by background effective stress**, while the **diurnal vs.
  semidiurnal variation is controlled by frictional properties and nucleation
  time.**

**This substantially validates hypothesis 9** — β does track effective normal
stress, as published work confirms. It also means the core idea is not ours to
claim.

*Strategic consequence.* Novelty must move to where the existing work is not:

1. **Long-period constituents.** Published work is diurnal/semidiurnal. Mf, Mm,
   Ssa, Sa, and the 18.61 yr nodal term are largely unexploited.
2. **Global scale.** Existing studies are single-fault (San Andreas, Cascadia). A
   global β field is genuinely new.
3. **The multi-basis / generalised Doodson extension.** Still unclaimed.
4. **Other bodies.** Heimisson & Avouac explicitly frame their model as applicable
   to "other solid-surface bodies" — supporting the moonquake plan.

This is a better position than it sounds: the mechanism is now *established
literature we can cite* rather than something we must prove.

### Finding 3 — pulsar timing array methods transfer cleanly ✓

Hypothesis 10 confirmed as viable. PTA analysis solves the same statistical
problem shape — small periodic signal, strongly red noise, rigorous false-alarm
control:

- Red noise modelled as a **Gaussian process in the Fourier domain** with
  harmonically related sinusoids
- **Analytic marginalisation** over individual Fourier component amplitudes while
  fitting the power-law amplitude and spectral index of the underlying process
- **Hierarchical Bayesian** treatment with hyperparameter marginalisation
- False-alarm null estimated **empirically**, because correlations make analytic
  nulls unreliable

Mature tooling exists (`enterprise`, and the EPTA/NANOGrav/MeerKAT analysis
stacks).

**Direct mapping to our problem:** catalogue red noise (aftershock clustering, Mc
drift, hydrological loading) becomes a fitted power-law Gaussian process; tidal
harmonics are the deterministic signal on top. Strictly more rigorous than the
Schuster test, which assumes independent samples and has no red-noise treatment.

Pleasingly, PTA's empirical false-alarm estimation is the same instinct as our
time-shifted null (doc 04 §6a) — independent convergence on the same safeguard.

### Finding 4 — the hydrological confounder is also the control

*Tidal and hydrological seismicity modulations reveal pore fluid diffusion during
earthquake nucleation* (**Science Advances**, sciadv.ady6350) treats tidal and
seasonal/hydrological modulation **jointly** to infer pore fluid diffusion.

Reframes doc 07's "most dangerous confounder": hydrological loading is a
**second periodic forcing with independently known amplitude and phase**. Rather
than merely controlling for it, use it as a **calibration signal** — it probes the
same transfer function at annual period, where tidal amplitude is weak.

Two forcings at different periods constrain the transfer function far better than
one. The confounder becomes an instrument.

### Not obtained this pass

- Ader et al. 2014 full text (cert failure on the Caltech mirror; paywalled at OUP)
- Heimisson & Avouac 2020 full text
- Science Advances paper full text (403)
- CTE / HW95 constituent catalogues

All four remain on the critical path. Acquiring the actual equations is the
blocking task for the next pass.

---

## 2026-08-21 — Third pass: primary sources obtained

**Access note.** No institutional access. Both critical-path papers were obtained
free from the author's public page at `web.gps.caltech.edu/~avouac/publications/`
(self-archived copies). **Author pages are the highest-yield OSINT route** — see
[11-osint-access.md](11-osint-access.md).

### The actual equations

**Heimisson & Avouac (2020), GRL 47** — general oscillatory stress:

```
Short period (T ≪ t_a):     R(t) = R₀ · exp( S_T(t) / (Aσ₀) )        (eq. 4, 8)
                            R₀ = r / M,   M = ⟨exp(S_T/Aσ₀)⟩ ≥ 1     (eq. 5, 7)

Long period  (T ≫ t_a):     R(t)/r ≈ 1 / (1 − t_a Ṡ_T(t)/(Aσ₀))     (eq. 9)

Characteristic time:        t_a = Aσ₀ / ṡ₀
```

with `S(t) = τ(t) − μσ(t)` the *modified* Coulomb stress, and μ = τ₀/σ₀ − α where
α is the Linker–Dieterich constant (0–0.25). **μ here is not a friction coefficient
in the traditional sense** — worth noting, since it is easy to misread.

**Ader et al. (2014), GJI 198** — finite rate-and-state fault:

```
Critical period:            T_a = 2π t_a = 2π aσ̄ / τ̇_a               (eq. 8)
Small-perturbation amplitude: R̃/r = Δτ / (aσ̄)                        (eq. B7)
```

### Finding 1 — band-pass CONFIRMED, and it comes from fault finiteness ✓

The second-pass reconstruction was right. Ader et al. state explicitly there is
"a critical period T_a, **at which the amplitude of the seismicity response
peaks**."

Structure of the finite-fault response:

| Regime | Amplitude | Phase |
|---|---|---|
| T < T_a | **increases with period** | ≈ 0, drifting to about **−π/4** as T → T_a |
| T ≥ T_a | falls; Coulomb-failure behaviour, rate ∝ stress rate | — |

The rising limb below T_a is a **finite-fault effect not present in the
spring–slider (SRM) model**, and so is the gradual −π/4 phase drift. Both are
therefore **discriminators between fault models**, exactly as hypothesis 2 wanted.

Also: finite-fault response amplitude is **always much larger than SRM predicts,
sometimes by more than an order of magnitude**. Consequence stated by the authors:
**aσ̄ values inferred from observations using SRM are systematically
underestimated.** Heimisson & Avouac note a second underestimation bias — Aσ₀ is
underestimated if stress is approximated by a single harmonic. Two independent
biases, both in the same direction.

### Finding 2 — hypothesis 9 was structurally wrong; the corrected version is better ✓

The synthesis claimed β ∝ τ̇/(Aσ₀). **That conflates two separate roles.**

```
Response AMPLITUDE      ←  Aσ₀ alone        R̃/r = Δτ/(aσ̄)
Peak/critical PERIOD    ←  T_a = 2π Aσ₀/τ̇   (sets which periods are favoured)
```

Stressing rate τ̇ does **not** set sensitivity. It only sets where the spectral peak
sits.

**This is better than the original claim**, because the two observables are
separable:

> Measure response **amplitude** → recover **Aσ₀** (effective normal stress × a).
> Measure the **peak period** T_a → recover **Aσ₀/τ̇**, hence **τ̇** independently.

**The harmonic spectrum yields effective normal stress and stressing rate as two
independent measurements.** And τ̇ — accumulating stress — is precisely the quantity
earthquake forecasting wants.

This also matches EPSL (2025) exactly: modulation *amplitude* controlled by
background effective stress; *frequency dependence* controlled by frictional
properties and nucleation time. Two observables, two parameters, confirmed in
published data.

### Finding 3 — numbers, and a resurrection of the long-period ranking

`T_a = 2π Aσ₀/τ̇`. Plugging in literature values:

| Setting | Aσ₀ | τ̇ | T_a |
|---|---|---|---|
| Ordinary crust (aftershock studies) | 0.01–0.1 MPa | ~0.003 MPa/yr | **~20–200 yr** |
| Parkfield tremor (Thomas et al. 2012) | 6 × 10⁻⁴ MPa | ~0.003 MPa/yr | ~1.3 yr |
| Tremor during slow slip (τ̇ ~100× higher) | 6 × 10⁻⁴ MPa | elevated | **~days** |
| Ader et al. simulation (σ̄ = 5 MPa) | — | — | ≈ 0.03 yr (~11 d) |

Two consequences:

**1. Tremor's T_a sits inside the tidal band.** That explains strong, clearly
frequency-dependent tidal modulation of tremor — and the fortnightly modulation —
without needing any special pleading. The mechanism for doc 05's "tremor is the
high-SNR testbed" is now quantitative.

**2. For ordinary crust, T_a ~ decades — so every tidal constituent, out to the
18.61 yr nodal term, lies in the T < T_a regime where amplitude RISES with
period.** The long-period ranking (nodal > Sa > Ssa > Mm > Mf > diurnal >
semidiurnal) is **plausible again** — but as the *rising limb toward a peak at
T_a*, not as the 1/ω low-pass falloff of the first pass. Same ordering, completely
different mechanism, and now with a theoretical amplitude curve to fit rather than
a bare slope.

### Finding 4 — two incidental confirmations

- **Ader et al. use Schuster spectra** (their Fig. 9) as the analysis tool. Our
  generalised-Schuster approach (doc 02) is the field's own method.
- **Harmonic forcing measurably modifies the magnitude distribution** in their
  simulations — whole magnitude bands "completely disappeared" under perturbation.
  Independent theoretical support for the Ide et al. (2016) b-value target
  (doc 04 §4).

### Non-linearity caveat

Response is `exp(S_T/Aσ₀)`, not linear. Linearisation holds when |S_T|/Aσ₀ ≪ 1:

- Ordinary crust: tidal stress ~10⁻³–10⁻⁴ MPa, Aσ₀ ~ 0.01–0.1 MPa → ratio 10⁻²–10⁻¹.
  **Linear is fine.**
- Parkfield tremor: Aσ₀ = 6×10⁻⁴ MPa → ratio **0.2–2. Strongly non-linear.**

So the tremor testbed needs the full exponential form, and M > 1 must be carried.
Note also **⟨R⟩ = r exactly** (eq. 6): oscillatory stress does not change the mean
rate, only the timing. Tides redistribute *when*, they do not create events. This
belongs in doc 00's claim discipline.

### RustSPICE update

Pulled to `c7f180b` (tag `v0.1.0`). Substantially restructured — see revised
[10-rustspice-requirements.md](10-rustspice-requirements.md).

---

## 2026-08-21 — Fourth pass: USGS public-domain set and tidal catalogues

**Access.** USGS Publications Warehouse and the USGS Rock Physics Lab page
(`earthquake.usgs.gov/research/eqproc/rockphysics/pubs.php`) delivered the entire
Beeler set free — USGS-authored work is public domain. Confirms doc 11's
prediction that this was the highest-value untried route.

Obtained: Beeler & Lockner (2003), Beeler et al. (2013) LFE rheology, Beeler et al.
(2016) effective stress, Lockner & Beeler (1999) premonitory slip.

### Finding 1 — the event-count requirement was wrong by two orders of magnitude ✓✓

**Previously recorded: 10⁵–10⁶ events. Actual: ~10⁴.**

That figure came from a search summary in the first pass and was never verified
against the paper. The paper itself, equation 18:

```
N ≥ ln(P_rw) / ( Δτ_u / (2 a σ_n) )²
```

**N scales as the inverse square of normalised stress amplitude.** Worked examples
from the paper:

| Stress amplitude | Events needed |
|---|---|
| Δτ = 0.01 MPa (a = 0.0045, σ_n = 18 MPa/km, 5–15 km) | **6.2 × 10³ – 5.5 × 10⁴** |
| Δτ = 0.1 MPa (10× larger) | **60 – 541** |

And from the abstract: daily Earth tides require **">13,000 earthquakes to
detect."**

**Consequence — the project is far more feasible than doc 01 recorded.** USGS
ComCat holds millions of events; 13,000 is trivially available, and many *regional*
subsets clear the bar comfortably. We are not power-limited at the global scale,
and the doc 04 §7 power analysis should be re-run with the correct law. The
"hard requirement" framing for the full catalogue pull stands, but for coverage and
confounder control rather than raw statistical power.

**Doc 01 corrected.**

### Finding 2 — the second timescale, and a startling band prediction ★★

Beeler & Lockner identify the **nucleation duration t_n**, distinct from Dieterich's
t_a. Their result:

> Fault strength's second-order dependence on sliding rate "determines the duration
> of nucleation and **damps the response to stress change at frequencies greater
> than 1/t_n**."

And, critically:

> "The experiments suggest that the **minimum typical duration of earthquake
> nucleation on the San Andreas fault system is 1 year.**"

This is the missing short-period cutoff — the second timescale the band-pass
reconstruction required (third pass). Combining with Ader's critical period:

```
responsive band:    t_n  <  T  <  T_a
                   ~1 yr  <  T  <  ~20–200 yr        (ordinary crust)
```

**Everything faster than about a year is damped by nucleation. Everything slower
than T_a rolls off.**

If this holds, the implications are severe and specific:

| Constituent | Period | Predicted |
|---|---|---|
| Semidiurnal (M2, S2, N2) | 12–12.4 h | **damped** |
| Diurnal (O1, K1, P1) | 24–26 h | **damped** |
| Fortnightly, monthly (Mf, Mm) | 13.7–27.6 d | **damped** |
| Semiannual (Ssa) | 182.6 d | marginal |
| **Annual (Sa)** | 365.26 d | **at the band edge** |
| **Lunar nodal** | 18.61 yr | **inside the band** |
| **LOD / decadal** | decades | **inside the band** |

**This is the strongest theoretical result the project has found.** It says the
responsive band for ordinary earthquakes is *years to decades* — orbital and
long-period timescales — and that the short-period tidal constituents everyone
studies are the ones physics predicts should NOT work.

It independently explains: why semidiurnal tidal correlation is weak everywhere;
why annual/monsoon signals correlate (Ader's Himalayan motivation); why Bendick &
Bilham's decadal LOD correlation is not absurd; and why tremor — with t_n and T_a
both far shorter — responds at tidal periods while earthquakes do not.

**It also justifies the project's long-period focus on physical grounds.** Not
because long cycles are meaningful, but because the fault itself low-passes
everything faster than a year.

⚠ **Caveat.** t_n ≥ 1 yr is an *extrapolation from laboratory experiments to the
San Andreas*, not a field measurement. It is the paper's own inference and should
be carried as a hypothesis with a large uncertainty, not a constant. But it is
published, mechanistic, and testable — precisely by measuring the response spectrum
(hypothesis 12).

### Finding 3 — HW95 already contains the planetary expansion ✗

**Doc 02's central novelty claim was wrong.**

The **HW95 catalogue** (Hartmann & Wenzel 1995, *GRL* 22) is a harmonic development
of the tide-generating potential containing **12,935 waves, of which 1,483 are due
to direct planetary effects** — Venus, Jupiter, Mars, Mercury, Saturn. Based on
DE200, spanning 1850–2150. Scope: Moon to degree 6, Sun to degree 3, planets to
degree 2.

So the "generalised Doodson expansion including planetary arguments" that doc 02
proposed as our novel contribution **was published in 1995.**

**This is a gift, not a setback.** We do not have to derive it — 1,483
rigorously-computed planetary tidal waves with correct amplitudes, produced by
geodesists, are exactly the multi-body harmonic feature basis the project wants.
Using an authoritative published basis is a far stronger position than deriving our
own.

**Corrected novelty claim.** Not the expansion. Rather:
1. Using HW95/KSM03 planetary waves as an **ML feature basis tested against
   seismicity** — not done
2. **Broadband response spectroscopy** across six decades (hypothesis 12)
3. A **global β field** — published work is single-fault
4. The **moonquake → tremor → earthquake validation ladder**

**Doc 02 corrected.**

### Finding 4 — the modern successor catalogue, and it is in our band

**Kudryavtsev & Cionco (2025)**, *New and updated long-periodic terms in harmonic
development of the Earth tide-generating potential*, arXiv:2508.18111, accepted in
*Geodesy and Geodynamics*. **Free on arXiv.**

- **38 terms with period longer than ~18 years**, amplitude ≥ 10⁻⁸ m²s⁻²
- A relatively large term at **~7.4 kyr**
- **Several NEW waves near the 18.61-year lunar nodal period**
- Built on **DE441**, planetary terms included, 30,000+ year spectral analysis
- Updates the **KSM03** expansion (Kudryavtsev 2004)
- Released in both **HW95 and KSM03 format**

This lands exactly in the band Finding 2 predicts is responsive. New waves near
18.61 yr are directly relevant. Highest-priority acquisition.

### Finding 5 — LFE effective stress is kilopascal-scale

Beeler et al. (2013), *Inferring fault rheology from low-frequency earthquakes on
the San Andreas*: observed tidal modulation "restricts ambient stress to be at most
a **few kilopascal**."

A few kPa = 10⁻³ MPa, consistent with Thomas et al.'s Aσ₀ = 6×10⁻⁴ MPa for
Parkfield tremor. Confirms the low-effective-stress mechanism for LFE/tremor tidal
sensitivity and supports tremor as the high-SNR testbed.

### Catalogue access status

**ETERNA 3.30** bundles seven tidal potential catalogues — Doodson (1921),
Cartwright et al. (1971, 1973), Büllesfeld (1985), Tamura (1987), Xi (1989),
Roosbeek (1996), Hartmann & Wenzel (1995). HW95 gives better than 1 nGal accuracy.

⚠ ETERNA 3.30 was distributed on CD-ROM for US$300 in 1997. **Modern free access
route not yet established.** Options to try next pass: ETERNA 3.4 / ETERNA-x
redistributions, IGETS / ICET (International Center for Earth Tides), PyETERNA and
other open reimplementations, and the KSM03/HW95-format catalogue released with
Kudryavtsev & Cionco (2025) — which may be the cleanest free path.

**RATGP95** (Roosbeek 1996) has a free PDF at Oxford Academic.

### Open

- ETERNA / HW95 catalogue data files — free source not yet found
- Kudryavtsev & Cionco (2025) catalogue data
- Ader PhD thesis (Caltech THESIS, open) — likely fuller derivations

---

## 2026-08-21 — Phase 1 validation: PASSED ✓

First end-to-end run. `cargo run --release --example moonquake_validation`.

**Data.** Apollo PSE expanded event catalogue, `levent.1008weber.csv`, from the PDS
Geosciences Node — public domain, no credentials. Fetch script at
`scripts/fetch-apollo.sh`.

**Ingestion finding — use `T2`, not `T1`.** The catalogue carries two
classification columns: `T1` original, `T2` revised (Nakamura 2005). They disagree
for 1,471 events, most of them originally logged as meteoroid impacts (`M`) and
later reclassified as deep moonquakes (`A`). Counting `T1` finds **1,359** deep
moonquakes; `T2` finds **7,082** — the figure the literature reports. Parsing to
usable times yields **6,954** (128 rows have unparseable time fields and are
skipped). Nest A1 has **424** events.

**Result — five of six known periodicities recovered, all nests (N = 6,954):**

| Period | Literature | Recovered | Error | Power |
|---|---|---|---|---|
| Half-month | 13.60 d | **13.609 d** | +0.07% | 869.4 |
| Draconic month | 27.212 d | **27.19 d** | −0.09% | 81.2 |
| Anomalistic month | 27.5546 d | **27.567 d** | +0.05% | 235.5 |
| Synodic month | 29.5306 d | **29.58 d** | +0.16% | 83.5 |
| Solar perturbation | ~206 d | **206.19 d** | +0.09% | 111.7 |
| Draconic/anomalistic beat | ~2190 d | no peak | — | — |

**Nest A1 (N = 424):** the draconic month is the **strongest** peak (power 131.5),
matching the literature's characterisation of A1. The 6-year beat appears at
2083 d (−4.87%, power 43.8) — acceptable given the record spans only ~1.4 cycles.

**Independent check of the "known" periods.** Both long periods derive from lunar
month beats, computed here from first principles:

```
draconic/anomalistic beat      = 2188.6 d    (literature ~2190, the 6-year cycle)
(anomalistic/synodic beat) / 2 =  205.9 d    (literature ~206)
```

So the 206-day period is **half the perigee–syzygy ("full moon") cycle**, and the
6-year period is the draconic/anomalistic beat. Neither is folklore.

### Unexplained peak — flag for the null

All-nests shows the **second strongest peak of the whole spectrum at 1886.6 d
(power 625.0)**, which matches no known deep moonquake period.

Most likely a **catalogue artifact**: the record spans ~2,920 days, Apollo stations
came online between 1969 and 1972 and ran to 1977, so overall event rate varies on
multi-year timescales and the observational envelope has its own spectral
structure. A 901 d peak (power 128.6) is similarly suspect.

**This is exactly what the time-shifted null exists to catch**, and it is a useful
early demonstration that raw periodogram power is not evidence. Next step for
Phase 1 is to run the null against these peaks and confirm the known periods
survive while 1886 d does not.

Smaller peaks at 6.805 d (= 13.609/2) and, in A1, 9.11 d (≈ 27.19/3) are harmonics
of the recovered periods, as expected.

### Assessment

**The instrument works.** The pipeline recovers five known tidal periodicities from
a real catalogue to better than 0.2%, on data where tidal forcing is the dominant
mechanism. Per `docs/12-build-plan.md` this was the gate on everything downstream.

Not yet done in Phase 1: the time-shift null on these peaks, and the ephemeris path
— computing actual tidal phase from Earth–Moon–Sun geometry rather than folding on
trial periods. The latter needs lunar kernels and is what carries forward to
Phases 2 and 3.

---

## 2026-08-21 — A3/A4: tidal phase implemented, and the null kills a 10⁻⁸⁹ result

**Code:** `ph-core::phase`, `examples/moonquake_tidal_phase.rs`.

### What was built

Tidal phase by the standard construction (Tanaka et al. 2002 and successors):
sample the forcing finely, locate successive **maxima**, assign each event a phase
by interpolating between the maxima bracketing it. Peak times refined by parabolic
interpolation for sub-sample accuracy.

Forcing scalar: largest eigenvalue of the lunar tidal tensor from real Earth + Sun
geometry in `MOON_PA`. Because the Moon is tidally locked, the tensor's
*orientation* barely moves and its *magnitude* varies as GM/d³ — so maxima are
lunar perigees and phase 0 means perigee.

**Self-check passed:** mean interval between forcing maxima came out **27.539 d**
against the anomalistic month's 27.555 d, confirming the scalar tracks lunar
distance as intended.

### The result

6,954 deep moonquakes, all phased, none dropped.

```text
analytic Schuster, order 1:   D²/N = 204.9   p = 1.07e-89   phase −92.1°
                     order 2: D²/N =  39.6   p = 6.52e-18   phase −91.7°

time-shift null, 308 offsets: observed D² = 1.42e6
                              null max    = 1.13e7      empirical p = 0.699
```

**The analytic test says p ≈ 10⁻⁸⁹. The time-shift null says p = 0.70.**

Nearly a third of shifted realisations cluster at least as strongly as the true
alignment, and the strongest cluster eight times harder. The analytic p-value
overstates significance by roughly **88 orders of magnitude**.

### Why — the degeneracy is broader than documented

The second-pass caveat said the time-shift null is degenerate against a *single
exact frequency*. That was too narrow. **Quasi-periodicity in the forcing is not
sufficient.** What matters is whether the **catalogue** shares the forcing's period.

Deep moonquakes are themselves locked near the anomalistic month. So is the
forcing. A global shift then rotates the phase cluster to a different phase while
leaving its concentration intact — exactly the invariance, arriving through the
catalogue rather than through the forcing.

Consequence: for this dataset, **"do events cluster in tidal phase?" is not a
falsifiable question.** Both series are periodic at the same rate, so clustering is
guaranteed at *some* phase.

### What this does and does not mean

**It does not mean deep moonquakes are not tidally driven.** The literature is
unambiguous that they are, and our own periodogram recovers five tidal periods to
better than 0.21%.

**It means this particular test cannot demonstrate it** — and, more importantly,
that reporting the analytic p-value alone would have been a false positive of
spectacular size. The null did its job. This is the first time in the project that
the machinery has caught something the naive analysis would have got badly wrong,
and it caught it on the dataset we chose *because* we knew the answer.

**The falsifiable question is which phase**, and whether it is consistent across
independent sub-populations in a way a failure mechanism predicts. Which is
precisely why Weber, Bills & Johnson (2009) work **per moonquake nest** rather than
on the pooled catalogue — a design choice whose necessity is now measured rather
than inherited.

### Consequences for the plan

- **A4 does not pass as specified.** The pooled-catalogue null cannot resolve the
  1886 d and 901 d artifacts either, for the same reason. A4 needs restating as a
  per-nest test.
- **A5/A6 rise in priority.** Per-nest Coulomb projection against Weber's cluster
  constraints is now the only route to a falsifiable Phase 1 claim, not merely a
  second validation.
- **Carry this to Phase 3.** Earthquake catalogues have strong internal
  periodicity too — aftershock sequences, seasonal detection cycles. The same trap
  is waiting, and the analytic Schuster p-value will be just as wrong there.
- `nakamura_2005_dm_locations.csv` (106 nests, already downloaded) becomes the
  next thing to use.

---

## 2026-08-22 — A5/A6: per-nest Coulomb search, and a second null artifact

**Code:** `ph-core::fault`, `examples/moonquake_nest_coulomb.rs`. 43 tests.

### What was built

Fault-plane resolution in the Aki & Richards convention: normal and slip vectors
from strike/dip/rake, traction decomposition, `ΔCFS = τ + μ′σₙ`, and rotation of a
body-fixed tensor into a local North-East-Down frame. Verified against invariants —
unit and orthogonal basis vectors, `|t|² = σₙ² + |τ|²`, rotation-invariant trace and
eigenvalues.

**Two things worth recording.**

*Depth does not matter.* The degree-2 tide-generating potential goes as
`r²·P₂(cos ψ)`, so its second derivatives are **constant throughout the body**. Deep
moonquake nests at 700–1200 km see the same tidal tensor as the surface. Only
latitude and longitude enter, and only to set the local frame. (Elastic response
does vary with depth, but that is the Love-number scale factor, which does not
change timing.)

*Coulomb stress is linear in the tensor*, `ΔCFS = T_ij C_ij` with
`C_ij = ½(uᵢnⱼ + uⱼnᵢ) + μ nᵢnⱼ`. So the variance of ΔCFS over any time series is a
quadratic form `cᵀΣc` on the **covariance of the tensor components**, precomputed
once. A grid search over thousands of planes and thousands of epochs collapses from
a full pass per plane to ~72 flops per plane — exactly, not approximately. This is
what made a 2000-trial null over 74 nests run in a minute.

### The result — and a null that lied

Weber's criterion: search fault orientations for the one minimising
`score = std(ΔCFS at events) / std(ΔCFS over the span)`. Search space 1,728 planes ×
5 friction values. Nests with ≥20 events: **74**.

| Null construction | Nests below p = 0.05 |
|---|---|
| **Uniform** random event times | **73 / 74** |
| **Shift** — slide the sequence, preserve spacing | **17 / 74** |

**The uniform null was wrong and would have produced a spectacular false
positive.** Drawing random times destroys the catalogue's temporal clustering, so
the test silently became "are events clustered in time at all?" — trivially yes,
since events at nearby epochs see nearly identical tensors and therefore low
variance, with no tidal alignment required.

The shift null preserves relative spacing and breaks only the alignment with the
forcing, which is the question actually intended.

**This is the second time a null has manufactured significance in this project**,
after the pooled-phase test (A3/A4). Both had the same root cause: *a null that
fails to preserve the catalogue's own structure tests a different hypothesis than
the one stated.*

### Honest final numbers

```text
shift null, 2000 trials:   17/74 nests nominally p < 0.05   (3.7 expected)
Benjamini-Hochberg, FDR 0.05:   0/74 nests survive
ensemble excess:                17 vs 3.7, sd 1.87  ->  7.1 sigma
```

**No individual nest is defensible after multiple-comparison correction. The
population-level excess is 7.1σ.**

That is a real but modest result, and it matches Weber, Bills & Johnson (2009)
directly: they report that for some clusters the constant-stress fit is good while
for others it is "not strongly dependent on plane orientation." A heterogeneous
population with a minority of well-constrained nests is precisely what 17/74 looks
like.

Note also that the p-value floor bit in practice. BH over 74 tests needs
p ≤ 6.8×10⁻⁴ for the strongest nest to survive; at 200 trials the floor was
5×10⁻³, so **no nest could have passed regardless of signal strength.** Raising to
2000 trials moved the floor to 5×10⁻⁴, below the threshold — after which the answer
was still zero, but now for a real reason rather than an arithmetic one.

### Caveats

1. The shift null may itself be partially degenerate, since catalogue and forcing
   share the anomalistic month (A3/A4). 17/74 could still be inflated.
2. No Love numbers, so this is stress **shape**, not magnitude. Timing and the
   orientation search are unaffected; absolute ΔCFS in Pa is not available and must
   not be reported from this code.
3. Grid resolution is 15°/15°/30°. A finer grid would lower every score, observed
   and null alike, so the comparison holds — but per-nest orientations should not
   be quoted at better than grid resolution.

### Plan status

**A5 complete. A6 complete, with a negative headline and a positive ensemble.**
Phase 1 has now produced its real deliverable, which is not "tides drive
moonquakes" — that was known — but a **validated instrument plus two documented
ways to fool it**. Both traps are waiting in Phase 3, where earthquake catalogues
carry aftershock clustering and seasonal detection cycles of their own.

---

## 2026-08-22 — C1/D1: Parkfield LFEs, and the M2/S2 gate fires

**Code:** `ph-core::parkfield`, `examples/parkfield_constituents.rs`. 48 tests.

### The dataset

Shelly (2017, updated 2024), from USGS ScienceBase — **public domain**.
**1,528,117 LFEs, 88 families, 23.1 years.** Every family holds ≥5,333 events;
median 16,058, max 44,156. Compare the largest Apollo moonquake nest at 85.

LFE **families** are the direct analogue of moonquake **nests**. Most families
individually clear Beeler & Lockner's ~10⁴ event requirement.

### The measurement

Schuster power evaluated at the **exact** periods of named constituents rather
than by blind peak search, which turns doc 08 §13b's validity gate into a direct
measurement. Pooled catalogue:

| Constituent | Period (d) | Power | vs S1 |
|---|---|---|---|
| M2 principal lunar | 0.5175 | 1,831 | 0.11 |
| **S2 principal solar** | 0.5000 | **7,729** | 0.48 |
| N2 larger elliptic | 0.5274 | 21 | 0.00 |
| **K1 luni-solar diurnal** | 0.9973 | **18,861** | **1.16** |
| **S1 solar diurnal** | 1.0000 | **16,245** | 1.00 |
| O1 principal lunar diurnal | 1.0758 | 8,052 | 0.50 |
| Mf lunar fortnightly | 13.661 | 334 | 0.02 |
| Msf lunar synodic fortnightly | 14.765 | 663 | 0.04 |
| Mm lunar monthly | 27.555 | 631 | 0.04 |
| Ssa solar semiannual | 182.62 | 431 | 0.03 |
| Sa solar annual | 365.26 | 4,654 | 0.29 |

Null expectation for `D²/N` is **1**. The largest family alone reproduces every
ratio to within a few percent.

### The gate fires

**S1 is exactly 24.000 h and has essentially no body-tide amplitude.** Its power of
16,245 is therefore a *direct measurement of the detection artifact floor* — LFE
detection is template matching on continuous data, and its sensitivity tracks the
day-night cultural noise cycle.

Three readings follow immediately:

1. **K1 (23.93 h) at 1.16× S1 is unusable**, exactly as doc 08 §13b predicted.
   Degenerate with the diurnal artifact.
2. **S2 (7,729) exceeds M2 (1,831) by 4.2×.** For a genuine body tide the ordering
   must be the other way — the M2 tidal potential is roughly 2.2× S2's. Seeing S2
   dominate is the signature of a **thermal or cultural semidiurnal cycle**, not a
   body tide. **The gate returns "suspect artifact."**
3. **Sa (4,654) is large but confounded** with seasonal variation in detection.

### What this does and does not mean

**It does not contradict the literature.** Thomas et al. (2012) find strong tidal
modulation of Parkfield LFEs — but they resolve *tidal stress on the fault* and
test its phase. They do not fold event times on trial periods.

**The lesson is sharper than "the catalogue is noisy":**

> **Raw period folding on a detection-limited catalogue measures the detector, not
> the Earth.**

Our own doc 02 argued this from theory. It is now measured, on a catalogue chosen
because the physical effect there is known to be strong.

**This is the third trap the machinery has caught**, after the pooled-phase null
(A3/A4) and the uniform-time null (A5/A6). All three share a shape: *a statistic
that looks decisive while silently answering a different question.*

### Consequence for the band prediction

**This run does not test it.** The long-period constituents (Mf 334, Msf 663,
Mm 631, Ssa 431) are all far below the diurnal artifacts — but that comparison is
meaningless while the diurnal band is artifact-dominated. Testing the 1 yr–200 yr
band prediction requires the artifact removed first.

### Next (C2)

Do it properly, as Thomas et al. do: compute the tidal tensor at the Parkfield
family locations, resolve ΔCFS onto San Andreas geometry (`ph-core::fault` already
does this), and test phase against a structure-preserving null. The full
exponential response `R₀exp(S_T/Aσ₀)` with `M > 1` is required here — Thomas et al.
infer `Aσ₀ = 6×10⁻⁴ MPa`, putting Parkfield firmly in the non-linear regime.

Also worth adding: an explicit **alias analysis** (doc 08 §13e). With a 24 h
detection cycle this strong, beats against it will appear elsewhere in the
spectrum and need blacklisting.

---

## 2026-08-22 — C2: first result that survives correction

**Code:** `examples/parkfield_coulomb_phase.rs`.

### The null, fixed before the run

Every trap so far came from choosing a null after seeing the data, so this one was
specified first and justified on structure.

D1 established that the catalogue carries a large detection artifact locked to
**solar time** (S1 power 16,245 against a null expectation of 1). So:

> **Null: shift every event time by a whole number of solar days.**

- The artifact is locked to local solar time → a whole-day shift leaves it
  **exactly invariant**.
- The lunar tide precesses ~50 min per solar day → a whole-day shift **does**
  decorrelate ΔCFS phase.

It holds the confound fixed while sliding only the quantity under test. 400 shifts
drawn from ±(30–4000) days, excluding near-multiples of the 27.55 d anomalistic and
29.53 d synodic months.

### Setup

ΔCFS from the real Earth tidal tensor (Moon, Sun, planets) in `IAU_EARTH`, rotated
to local NED at 35.635 N, −120.150 E, resolved onto deep San Andreas geometry —
strike 137°, dip 90°, rake 180°, μ = 0.4. Sampled at 0.02 d over 23 years plus
padding: **841,272 points**.

**Self-check passed:** mean interval between ΔCFS maxima came out **0.5175 d** —
exactly M2. The forcing is what we think it is.

### Result

| | |
|---|---|
| Families tested (largest 12) | 12 |
| Nominally p < 0.05 | **10 / 12** (0.6 expected) |
| **Surviving Benjamini-Hochberg, FDR 0.05** | **9 / 12** |
| `D²/N` range | 106 – 1702 (null expectation 1) |
| Preferred phases | **all within 71.4°**, spanning 0.7° to 72.1° |

**The null has demonstrable power.** Two families are *not* significant, with
p = 0.14 and p = 0.22, and observed `D²` sits near — not far above — the null
maximum for most families. That spread is the signature of a working null,
in direct contrast to A5/A6's 73/74 pinned at the floor.

Every family clusters **shortly after the ΔCFS maximum**, consistent with Thomas
et al. (2012) for these same LFEs.

Geometry is robust: `D²/N` moves only 404 → 415 across strike ±10° and dip 80–90°.
A thrust geometry gives both lower power (356) and a very different phase
(−108.8°), so the statistic does discriminate between orientations.

### Caveats — none of these are small

1. **The 12 families are co-located within ~30 km**, so they see nearly identical
   forcing. Phase coherence across them is a **coherence check, not independent
   confirmation**. It shows the pipeline returns a consistent answer; it is not 12
   independent tests of phase.
2. **Absolute phase is geometry-dependent.** Strike 127 → 137 → 147 shifts the
   preferred phase 57° → 41° → 25°. The *existence* of clustering is robust; the
   *value* should not be quoted without stating the assumed geometry.
3. **Effect sizes are modest.** For most families observed `D²` is at roughly the
   98th–99th percentile of the null, not orders above it.
4. **No Love numbers**, so this is stress *shape*. Timing is unaffected; absolute
   ΔCFS in Pa is unavailable.
5. **This does not test the band prediction.** M2 is a 12.42 h constituent, deep
   inside the range doc 07 predicts should be *damped* for ordinary crust. That it
   works at Parkfield is consistent with tremor's `T_a` sitting in the tidal band —
   but the 1 yr–200 yr prediction remains **entirely untested**.

### Where this leaves the project

Phase 1 delivered a validated instrument and three documented traps. **C2 is the
first time the instrument has been pointed at contested data with a
pre-specified null and returned a result that survives correction.**

That is the methodological milestone the validation ladder existed to reach: not
"we found a correlation," but "we found one *after* three separate opportunities to
fool ourselves, using a null fixed in advance."

---

## 2026-08-22 — C3: the amplitude law, and a fourth trap

**Code:** `ph_core::phase::cycle_amplitude_at`, `examples/parkfield_amplitude_law.rs`.
50 tests.

### Why amplitude rather than constituent

D1 showed the diurnal and semidiurnal bands are dominated by a detection artifact
locked to solar time, so a constituent-by-constituent spectrum runs straight back
into it.

But **the artifact does not care how strong the tide is.** Binning events by the
peak-to-trough amplitude of the ΔCFS cycle they fall in gives a test that detection
bias cannot produce, because it is amplitude-independent by construction.

This is a direct test of the amplitude law `R̃/r = Δτ/(aσ̄)` (Ader et al. 2014, eq.
B7). Predictions fixed before the run: flat → artifact; `D²/N ∝ amplitude²` →
linear response; steeper → the non-linear exponential regime.

### Result

All 1,528,117 events, six equal-count amplitude bins:

| Bin | Events | Amplitude (rel.) | `D²/N` | ratio to bin 1 |
|---|---|---|---|---|
| 1 | 254,686 | 1.00 | 216 | 1.00 |
| 2 | 254,686 | 1.48 | 2,301 | 10.67 |
| 3 | 254,686 | 1.85 | 3,513 | 16.30 |
| 4 | 254,686 | 2.20 | 3,883 | 18.01 |
| 5 | 254,686 | 2.53 | 7,639 | 35.44 |
| 6 | 254,687 | 3.00 | 15,288 | 70.91 |

**Monotonic across every bin. A 3× amplitude range produces a 71× change in phase
concentration.** Log-log slope **3.56**, against 2 for a linear response.

### The fourth trap — per-bin nulls are compromised

Bin 6 has the **highest** `D²/N` (15,288) and yet **p = 0.45**.

That is the tell. **Binning by amplitude selects on the forcing itself**, and
high-amplitude cycles recur at the ~14.77 d spring–neap period, so every bin
inherits temporal structure derived from the very signal under test. A whole-day
shift does not preserve that structure, so the per-bin null answers a different
question — again.

**The per-bin p-values in the table above are not trustworthy and should not be
quoted.**

### Nulling the claim actually being made

The claim is the *trend*, so the trend is what must be nulled: shift the event
times, **re-bin at the shifted times**, refit the slope.

```text
observed slope 3.56    null median 0.66    null max 3.53    p = 0.0050
```

**The trend survives.** Null slopes centre near 0.66 — no trend, as expected when
alignment is broken — while the observed 3.56 exceeds all 200 draws.

⚠ **But only just.** Null max 3.53 against observed 3.56 is a thin margin, and
p = 0.0050 is the 200-trial floor rather than a measurement. The true value is
plausibly somewhat higher. More trials would sharpen it; the null slope
distribution has a heavy tail and deserves more sampling before this is published.

### Reading

Response scales with tidal forcing amplitude, **faster than linearly**. That is
consistent with the exponential form `R = R₀exp(S_T/Aσ₀)` being in its non-linear
regime — which is exactly what Thomas et al. (2012) imply for Parkfield with
`Aσ₀ = 6×10⁻⁴ MPa`, giving `S_T/Aσ₀ ≈ 0.2–2`.

Two independent lines now agree that Parkfield LFEs respond to tidal Coulomb
stress: C2's phase clustering (9/12 families surviving FDR) and C3's amplitude
scaling. They rest on different statistics and different nulls.

### Trap count: four

1. Pooled-phase null degenerate when catalogue shares the forcing's period
2. Uniform-time null tests temporal clustering, not tidal alignment
3. Raw period folding measures the detector
4. **Per-bin nulls compromised when the binning variable derives from the forcing**

All four share one shape: *a statistic that looks decisive while silently answering
a different question.* Each was caught only because the result was checked against
what the null was actually testing, rather than against intuition.

### Still untested

**The band prediction.** Every result so far sits at M2 (12.42 h) — deep inside the
range doc 07 predicts should be *damped* for ordinary crust. Parkfield works
because tremor's `T_a` sits in the tidal band. The 1 yr–200 yr claim for ordinary
crust remains untouched, and needs a terrestrial earthquake catalogue (Phase 3).

---

## 2026-08-22 — C3b confirms; C4 fails, and the failure is instructive

### C3b — the amplitude law holds at 2,000 trials

```text
observed slope 3.56   null median 0.43   null max 4.09   p = 0.0095
```

Against the 200-trial run (p = 0.0050, null max 3.53). As flagged, **that p was a
floor, not a measurement** — with ten times the sampling, 18 of 2,000 null draws
now exceed the observed slope and the null maximum overtakes it.

**The result survives at p ≈ 0.01** rather than 0.005. Weaker than it first looked,
still significant, and now properly estimated. The prediction that the floor was
hiding the true value was correct.

Also required an optimisation worth noting: sorting 1.5M events per trial would
have dominated runtime, so amplitude bin edges are now fixed once from the observed
distribution and reused. Binning is O(n) and — more importantly — **identical
between observed and null**, which is what a like-for-like comparison needs.

### C4 — frequency-resolved response: method failed

Built `ph_core::demod` (complex demodulation, 5 tests, verified to separate M2 from
a 4× larger S2 when the window exceeds their 14.77 d beat). Then attempted response
per constituent at Parkfield.

**Two nulls tried, both invalid.**

**Trap 5 — the time-shift null, recurring.** Demodulation isolates one constituent,
which makes the band a **near-pure tone** — and against a pure tone a time shift
merely *rotates* the phase cluster without diluting it. `D²` is near-invariant and
the null has no power.

This is trap 1 exactly. `ph_core::stats` documents it in a module-level warning.
**The demodulation step recreated the degenerate case, and I walked into it
anyway.** Worth recording plainly: documenting a trap is not the same as being
immune to it, because the trap reappears wearing different clothes.

**Trap 6 — the sham-frequency null.** Replacing it with "run the same procedure at
frequencies carrying no tide" seemed sound. The result:

```text
sham floor (10 quiet bands):  median D²/N = 1,439,329   max 1,520,587
M2:  17,986      O1: 8,357      Mf: 4,195      Sa: 1,488,398
```

With N = 1,528,117, a sham median of 1.44M means **D² ≈ N²** — events land at
essentially *one* phase at frequencies where there is no tide. Every real
constituent sits **below** the floor.

The diagnosis: at a frequency with no genuine power, `z̄` is dominated by whatever
leaks in from the nearest strong constituent at ω′. Then `z̄ ∝ e^{i(ω′−ω)t}`, so
`arg z̄ = (ω′−ω)t`, and the reported band phase `ωt + arg z̄` collapses to `ω′t` —
**the phase of the leaking constituent, not of the target.** A tide-free frequency
does not give a neutral baseline; it gives a relabelled copy of the dominant band.

### What C4 needs instead

**Do not derive constituent phase by demodulating the composite ΔCFS series.**
Compute it from the **astronomical argument** — the exact Doodson combination of
fundamental arguments for that constituent. That phase is uniform over long spans
by construction, so events uniform in time are uniform in phase, and both the
statistic and its null behave.

This is what the tidal literature does, and it is why the field works in terms of
constituent arguments rather than filtered series.

**It also makes D2 a blocker.** Analytic constituent arguments mean the HW95/KSM03
expansion, which doc 07's fourth pass left as "free access not yet established."
D2 moves from parallel task to **prerequisite for C4**.

### Status

- **C4 is not done.** No frequency-resolved response exists yet, and no transfer
  function. The numbers above are artifacts and must not be quoted.
- `ph_core::demod` is correct and tested for what it does — isolating a constituent
  from a series. The error was using its *phase output* as a statistic against a
  null, not the demodulation itself. It stays.
- **C2 and C3 are unaffected.** Both used the full composite ΔCFS, which is
  genuinely quasi-periodic, so their nulls retain power. C3b independently
  confirmed at 2,000 trials.

### Trap count: six

1. Pooled-phase null degenerate when catalogue shares the forcing's period
2. Uniform-time null tests temporal clustering, not tidal alignment
3. Raw period folding measures the detector
4. Per-bin nulls compromised when the binning variable derives from the forcing
5. **Time-shift null degenerate against a demodulated single constituent** — trap 1
   recurring, walked into despite being documented
6. **Sham-frequency null invalid**: demodulating where there is no tide returns the
   leaking constituent's phase, not a neutral baseline

Six for six, all the same shape: *a statistic that looks decisive while silently
answering a different question.*

---

## 2026-08-22 — D2/C4 redone: the first transfer function

### D2 turned out not to need a download

The blocker was thought to be catalogue access. It was not. Constituent **phase**
is an integer combination of six astronomical arguments, each an analytic
polynomial in time — no catalogue required. Only constituent **amplitudes** need
HW95/KSM03, and phase is what the statistics use.

`ph_core::doodson` implements the six fundamental arguments (τ, s, h, p, N′, pₛ)
and 13 named constituents. **All six tests passed first run**, including the one
that matters:

- Every constituent period matches its published value to better than 0.1%,
  measured from the argument's rate rather than assumed
- **Phase histograms are flat to <5% over 12 bins** (M2) and <10% (long period) —
  *precisely the uniformity that demodulated phase lacked*
- S2 repeats exactly every half solar day; M2 does not

### The null that finally works

A **global** time shift can never work for a single constituent: `D²` is invariant
under rotation, and shifting one constituent's phase globally *is* a rotation. That
is a property of the statistic, not a fixable detail.

So the null shifts **each block independently**, with block length `max(4×period,
30 d)`. Within-block clustering is preserved; alignment between blocks is
randomised. Concentration changes, which is what a null must do.

### Result — 1,528,117 events, 23 years

| Band | Period (d) | Blocks | `D²/N` | Null median | Ratio | p |
|---|---|---|---|---|---|---|
| **M2** | 0.5175 | 280 | 17,975 | 86 | **208** | 0.0025 ✱ |
| **N2** | 0.5274 | 280 | 1,559 | 32 | **49** | 0.0025 ✱ |
| **O1** | 1.0758 | 280 | 8,020 | 65 | **123** | 0.0025 ✱ |
| **Q1** | 1.1195 | 280 | 254 | 31 | **8.3** | 0.0025 ✱ |
| Mf | 13.661 | 154 | 338 | 235 | 1.4 | 0.33 |
| Msf | 14.765 | 142 | 660 | 218 | 3.0 | 0.13 |
| Mm | 27.555 | 76 | 638 | 262 | 2.4 | 0.20 |
| Ssa | 182.62 | 11 | 431 | 666 | 0.6 | 0.65 |
| Sa | 365.26 | 5 | 4,654 | 1,783 | 2.6 | 0.12 |

**4/9 survive Benjamini-Hochberg at FDR 0.05 — and all four are short-period.**

**The response is band-limited, with the band at hours-to-a-day.** It falls roughly
two orders of magnitude between 0.5 d and 14 d, and every long-period constituent
is non-significant.

### Why this is the measurement we wanted

Doc 07's band prediction says **ordinary crust** should respond at **1 yr–200 yr**
and be damped below that. Parkfield is *tremor*, where `T_a` is short — so the
prediction for Parkfield is the **opposite** shape, and that is what we measure.

**Phase 3 now has a control.** If terrestrial earthquakes show the *same* shape as
Parkfield, the band prediction is wrong. If they show the mirror image — nothing at
M2, response at Sa and longer — it is confirmed. Either way the comparison is
decisive, which no single measurement could be.

### An internal check, half passing

`D²/N ∝ amplitude²`, so same-band constituent ratios test the amplitude law
independently of frequency (doc 08 §13a). Both pairs have amplitude ratios near
5.3, predicting a response ratio near 28:

- **O1/Q1 = 31.6** against 28 predicted — close
- **M2/N2 = 11.5** against 28 predicted — off by 2.4×

Partial support. Worth understanding rather than glossing: N2 responding more than
its amplitude warrants is the kind of discrepancy that either indicates a real
frequency dependence within the semidiurnal band or a leak we have not found.

### Caveats

1. **Long-period constituents are underpowered.** Ssa gets 11 blocks and Sa only 5.
   Non-significance there **cannot be distinguished from insufficient power**, and
   must not be read as evidence of absence.
2. The null median rises with period (86 → 1,783) because longer blocks mean fewer
   randomisations. Ratios are the fair comparison; raw `D²/N` across constituents
   is not.
3. No Love numbers, so this is response versus *frequency*, not response per unit
   *stress*. Converting to a true `R(ω)` still needs F2.

### Status

**C4 done.** `ph_core::demod` remains — correct for isolating a constituent, just
not for producing a null-valid phase. D2 is closed as a blocker.

Trap count stays at six. This run introduced none, which is the first time.

---

## 2026-08-22 — F2, and a correction to the previous entry

### F2 — elastic calibration

`ph_core::love` converts tidal tensors (s⁻²) to stress (Pa) via the standard
degree-2 surface relation: potential from the tensor's radial component, strain via
`(2h₂ − 6l₂)`, stress via Hooke's law. **Not** a full elastic solution — an
order-unity-accurate scalar calibration, documented as good to ~2×. 7 tests.

Validation, and a striking one:

```text
M2 solid Earth tide:  strain 9.92e-9    stress 595 Pa  (5.95e-4 MPa)
Thomas et al. (2012) Parkfield A*sigma0:      600 Pa  (6.0e-4 MPa)
```

**Independently calibrated M2 stress matches the published `Aσ₀` to 1%.** So
`S_T/Aσ₀ ≈ 1` at Parkfield — squarely in the 0.2–2 range Heimisson & Avouac
identify as non-linear. That is an independent explanation for **C3's slope of
3.56**, steeper than the linear prediction of 2, arrived at from elastic constants
rather than from the seismicity.

### ✗ Correction: C4 is not band-limited

The previous entry read C4's raw `D²/N` — 17,975 at M2 against 338 at Mf — as a
transfer function falling with period, and concluded the response was
"band-limited at hours-to-a-day."

**That was wrong.** `D²/N` is response to *whatever forcing exists at that
frequency*, and the tidal potential's own amplitudes fall with period. Dividing by
the ΔCFS amplitude at each constituent — obtained by least-squares regression on
the analytic Doodson phase, so it is the stress *this fault at this site* actually
sees — gives:

| Band | Period (d) | Amp (Pa) | `D²/N` | **R(ω) per Pa** | p |
|---|---|---|---|---|---|
| M2 | 0.5175 | 679 | 17,975 | **3.19e-4** | 0.0025 ✱ |
| N2 | 0.5274 | 130 | 1,559 | **4.92e-4** | 0.0025 ✱ |
| O1 | 1.0758 | 137 | 8,020 | **1.06e-3** | 0.0025 ✱ |
| Q1 | 1.1195 | 26 | 254 | **9.87e-4** | 0.0025 ✱ |
| Mf | 13.661 | 74 | 338 | 4.02e-4 | 0.33 |
| Mm | 27.555 | 37 | 638 | 1.12e-3 | 0.20 |
| Ssa | 182.62 | 34 | 431 | 9.93e-4 | 0.65 |
| Sa | 365.26 | 4.9 | 4,654 | 2.26e-2 | 0.12 [few blocks] |

**Response per unit stress varies by a factor of ~3 across 0.5 d to 27 d — a
50-fold range in period.** The 208× spread in raw ratio was almost entirely the
forcing amplitude spectrum.

Consequences:

1. **No band limitation is detected between 0.5 d and 27 d.** The response looks
   effectively scale-free there.
2. **`T_a` is not located.** The previous entry's `T_a ≲ 1 d` bound is withdrawn,
   and with it the "tension" between that bound and published `Aσ₀` — there was no
   tension, only a normalisation error.
3. **Non-significance at Mf, Msf and Mm is a power problem, not a physics one.**
   Their forcing amplitudes are small, so ε is small, so the same `R` yields less
   detectable concentration.
4. Sa's `R` of 2.26e-2 is 20× everything else, but it is non-significant with 5
   blocks and only 4.9 Pa of forcing. **`R` is unreliable wherever amplitude is
   small and response is non-significant**, since ε is then noise divided by a
   small number. It is not a result.

### The methodological point

A raw response measurement is not a transfer function. **Dividing by the forcing is
not a refinement, it is the definition** — and skipping it produced a confident,
wrong, physically-flavoured conclusion that survived one whole entry before the
normalisation caught it.

This did not come from a bad null. All the statistics were sound. It came from
reporting `D²/N` as though it answered a question it does not answer — which is the
same shape as the six traps, arriving through interpretation rather than through
the statistic.

### What this does to the plan

**Phase 3's control is stronger than expected.** If Parkfield's `R(ω)` really is
flat over 0.5–27 d, then a *flat* result for ordinary crust would say the band
prediction is wrong, while a rising `R(ω)` toward years would confirm it. The
comparison no longer depends on our reading of a roll-off we never measured.

**Still needed for the long-period end:** the non-significance at Mf/Mm/Ssa/Sa is
power-limited, and more events will not fix Sa's 5 blocks. A longer catalogue or a
second site is the only route.

---

## 2026-08-22 — C5: independent replication at Cascadia

**Code:** `ph_core::cascadia`, `scripts/fetch-cascadia.sh`,
`examples/two_site_comparison.rs`. 72 tests.

### Why this was the most important gap

Every Parkfield result rested on **one location with co-located families**. C2's
phase coherence across 12 families was explicitly labelled a coherence check rather
than confirmation, because they all see the same forcing. A shared instrumental
artifact remained a live explanation.

### The second site

**678,084 tremor detections, 2009–2024**, from the PNSN interactive tremor
catalogue (Wech). Public, no credentials.

| | Parkfield | Cascadia |
|---|---|---|
| Setting | strike-slip transform | subduction megathrust |
| Location | 35.6 N, 120.2 W | 40–50 N, 122–125 W |
| Events | 1,528,117 LFEs | 678,084 tremor |
| Span | 23.1 yr | 15.4 yr |
| **Detection** | **template matching** | **envelope cross-correlation** |

The detection difference is what gives the comparison force. Parkfield's diurnal
artifact (D1: S1 power 16,245 against a null expectation of 1) arises from template
matching against a time-varying noise floor. A different pipeline should carry a
*different* artifact.

### Result

Phase-only, so no fault geometry is involved and the sites are directly comparable.
Same analytic Doodson phases, same per-block shift null, 400 trials.

| Band | Period (d) | PK ratio | PK p | CS ratio | CS p |
|---|---|---|---|---|---|
| **M2** | 0.5175 | **223** | 0.0025 ✱ | **100** | 0.0025 ✱ |
| **N2** | 0.5274 | **45.7** | 0.0025 ✱ | **7.8** | 0.0050 ✱ |
| **O1** | 1.0758 | **137** | 0.0025 ✱ | **51** | 0.0025 ✱ |
| Q1 | 1.1195 | 6.6 | 0.0075 ✱ | 0.9 | 0.56 |
| Mf | 13.661 | 1.5 | 0.36 | 0.4 | 0.74 |
| Msf | 14.765 | 2.9 | 0.11 | 2.7 | 0.18 |
| Mm | 27.555 | 2.3 | 0.18 | 2.5 | 0.17 |
| Ssa | 182.62 | 0.6 | 0.72 | 1.7 | 0.28 |
| Sa | 365.26 | 2.8 | 0.087 | 2.2 | 0.11 |

**8/9 constituents give the same verdict. M2, N2 and O1 are significant at both.**

The single disagreement is **Q1** — the weakest constituent, significant at the site
with 2.3× more events and not at the other. That is where power fails first, and it
is the expected failure rather than a contradiction.

### What this buys

**The shared-artifact explanation is now very hard to sustain.** Template matching
and envelope cross-correlation would both have to manufacture the same
M2/N2/O1-significant, long-period-null pattern, at different latitudes, on
different fault geometries, over different epochs.

This is the first result in the project that is **replicated rather than
internally consistent**, and it retroactively strengthens C2, C3 and C4.

### An honest note on effect size

Cascadia's ratios run systematically lower (M2 100 vs 223; O1 51 vs 137). Part is
sample size — 678k against 1,528k — but not all of it. If ε were equal, `D²/N`
would scale with `N`, predicting a Cascadia/Parkfield ratio of 0.44; the observed
M2 ratio is 0.15, implying Cascadia's fractional modulation is roughly **0.6× that
of Parkfield.** Different setting, different sensitivity. Worth quantifying properly
once Cascadia has fault geometry and an `R(ω)` normalisation.

### Two API behaviours worth recording

`scripts/fetch-cascadia.sh` documents both, because each would silently corrupt a
catalogue:

1. **The API caps a response at 20,000 events and does not say so.** A yearly
   request returns exactly 20,000 with HTTP 200 and no truncation flag. The first
   fetch produced a tidy-looking 314,569 rows that were quietly wrong. Monthly
   chunking plus a near-cap warning fixes it.
2. **Empty windows return 404, not an empty result.** The catalogue starts mid-2009,
   so early months 404 legitimately and `curl --fail` aborts the run.

Neither is a trap in the statistical sense, but both belong in the same family:
*a result that looks complete while silently being something else.*

### Still open at the long-period end

Mf through Sa remain non-significant at both sites. Cascadia does not resolve this
— it has a **shorter** span (15.4 yr against 23.1), so it has *less* long-period
power, not more. Testing the band prediction still needs ordinary crust and a long
catalogue, which is Phase 3.

---

## 2026-08-22 — P3.4: the band prediction test. Inconclusive, with bounds.

**Code:** `ph_core::comcat`, `scripts/fetch-comcat.sh`,
`examples/band_prediction_test.rs`. 81 tests.

### Magnitude of completeness, settled empirically

Mc drift projects onto long-period features and manufactures exactly the signal
the band prediction looks for, so the threshold is the whole ballgame. Decadal
counts decide it:

| Threshold | 1970s | 1980s | 1990s | 2000s | 2010s | verdict |
|---|---|---|---|---|---|---|
| M5.0+ | 13,581 | 16,025 | 14,660 | 17,234 | 18,469 | **+36%, monotonic — incomplete** |
| **M5.5+** | 4,377 | 4,384 | 4,865 | 5,160 | 4,898 | **18% spread, no trend — stable** |
| M6.0+ | 1,127 | 1,287 | 1,535 | 1,585 | 1,494 | +33% |

**M ≥ 5.5, 1970–2025: 25,962 events over 55 years.** Above Beeler & Lockner's
~13,000 requirement.

Also added `phase_at_longitude`: tidal phase is *local*, and for semidiurnal
constituents a 180° longitude error is a full cycle. Without it a global catalogue
would be scattered across every phase and any real signal erased.

### Result

| Band | Period (d) | PK ratio | CS ratio | **EQ ratio** | **EQ p** |
|---|---|---|---|---|---|
| M2 | 0.518 | 251 ✱ | 92 ✱ | **0.6** | 0.68 |
| N2 | 0.527 | 50 ✱ | 8.1 ✱ | **2.7** | 0.17 |
| O1 | 1.076 | 116 ✱ | 47 ✱ | **0.2** | 0.88 |
| Q1 | 1.120 | 7.4 ✱ | 0.8 | 2.0 | 0.20 |
| Mf–Sa | 13.7–365 | 0.6–3.3 | 0.4–2.8 | 0.5–1.0 | 0.48–0.72 |

**Nothing is significant for ordinary crust at any constituent.**

### The bounds, which are the actual result

A null result is worth nothing without a detection limit. Converting to fractional
rate modulation `ε = 2√((D²/N)/N)`:

| Band | EQ upper bound | PK observed | Tremor exceeds bound by |
|---|---|---|---|
| **M2** | **< 3.88%** | 21.69% | **5.6×** |
| **O1** | **< 4.33%** | 14.49% | **3.3×** |
| **N2** | **< 4.14%** | 6.39% | **1.5×** |
| Mf–Sa | < 4.95–6.27% | 3.0–11.0% (ns) | — |

**Ordinary crust responds at least 3–5× less than tremor at M2 and O1.** That is a
real, quantified, defensible statement, and it is consistent with tremor's
anomalously short `T_a`.

### But the band prediction is still untested — say so plainly

The prediction is that ordinary crust responds **at Sa and longer**. Our
long-period bounds are only 5–6%, which is nowhere near tight enough to exclude a
modest long-period response. We have no positive long-period detection and cannot
rule one out.

**"Ordinary crust is quieter than tremor at short periods" is consistent with
damping, but it is equally consistent with ordinary crust simply having a larger
`Aσ₀` and no band structure at all.** The two are not distinguished by this data.

### What would settle it

To reach 1% sensitivity: `N = 4·D²/N_threshold / ε²` ≈ 4 × 10 / 10⁻⁴ =
**~400,000 events**. We have 25,962 — **the global M5.5+ catalogue is ~15× too
small.**

The trade is unavoidable. A regional catalogue at lower Mc (California M2.5+ since
1980 gives ~5×10⁵) reaches the count but sacrifices span and globality; the global
catalogue has the 55-year span the long-period end needs but not the events. There
may be no single catalogue that does both, in which case the band prediction is
answerable only by combining regions with per-region Mc control.

### One flagged coincidence, deliberately not claimed

**Sa is at p = 0.0973 at Parkfield and p = 0.0998 at Cascadia** — two independent
sites landing within 0.003 of each other at the annual constituent. Fisher's method
on the two gives p ≈ 0.055.

Not significant, and **Sa is precisely where seasonal hydrological loading
confounds** (doc 08 §11). Recorded because it is the only long-period feature
showing consistency across sites, and because pre-registering it now is the only
way to test it later without it becoming a post-hoc find. **It is not a result.**

### Status

**P3.4 attempted, inconclusive, bounded.** Everything upstream stands: the
moonquake validation, the two-site tremor replication, the amplitude law, the
elastic calibration. The central question remains open, and now has a number
attached to what would close it.

---

## 2026-08-22 — P3.5: depth stratification fails, and that is a stop signal

**Code:** `examples/depth_stratified_test.rs`.

### The pre-registered prediction

P3.4's global catalogue mixes every depth and fault geometry, diluting any effect
concentrated in one setting. Métivier et al. (2009) report the tidal anomaly is
**larger for shallower earthquakes**; Cochran et al. (2004) find their factor-3
effect in *shallow thrust* faults specifically.

So: **shallow responds more strongly than deep.** One split at 70 km, fixed before
running, no scanning over cut depths.

### Result — not supported

20,770 shallow (≤70 km) against 5,192 deep.

| Band | shallow ε | shallow limit | p | deep ε | deep limit | p |
|---|---|---|---|---|---|---|
| M2 | 0.13% | <4.80% | 0.99 | 4.18% | <6.92% | 0.14 |
| N2 | 2.07% | <4.22% | 0.17 | 1.52% | <7.46% | 0.78 |
| O1 | 1.04% | <5.16% | 0.72 | 1.60% | <7.98% | 0.77 |
| Mf | 2.52% | <5.57% | 0.39 | 2.17% | <7.53% | 0.59 |
| Mm | 1.14% | <6.48% | 0.90 | 2.44% | <7.99% | 0.51 |
| Sa | 3.13% | <7.49% | 0.59 | 3.47% | <8.35% | 0.39 |

**Nothing significant anywhere.** Shallow exceeded deep in **2/6 bands — fewer than
the 3 expected by chance.** Sign test p = 0.89 one-tailed. At M2 the shallow set is
essentially flat (0.13%) while the deep set is *higher* (4.18%), the opposite of
the prediction, though not significant.

### Why this is a stop signal, not a prompt to slice further

Splitting halved the sample, so every bound loosened — 4.8–8.4% here against
3.9–6.3% unsplit. **That was flagged in the example's own header before running,
and it is the general shape of the problem: each additional stratification makes
the bound worse, not better.**

Two stratifications in, nothing found, and the bounds are degrading. Continuing to
slice — by magnitude band, by region, by mechanism proxy, by cut depth — is exactly
the failure mode the pre-commitment in doc 16 §decision-4 exists to prevent:
**keep cutting until something crosses p < 0.05.** With six constituents and a
handful of strata, something eventually will.

**Stop slicing. The catalogue is power-limited and no partition fixes that.**

### What genuinely differs from slicing

One improvement is not a partition and is worth doing: **the earthquake analysis
used a weaker feature than it needed to.**

Both the tremor and earthquake tests used analytic Doodson phase — appropriate, and
apples-to-apples, since D² is unaffected by the constant longitude offset at a
single site. But raw tidal phase is blind to whether the tide actually *loads or
unloads each fault*. Parkfield's C2 result used ΔCFS resolved on a known fault
plane and found phase clustering; the global test has no mechanisms, so
compressional and extensional responses cancel in the pooled statistic.

**GCMT focal mechanisms would let us compute real ΔCFS per event** (`ph_core::fault`
already does the projection). That raises the signal by aligning the feature with
the physics rather than by discarding data — the opposite of stratification.

That is P3.2, and it is now the highest-value next step.

### Standing position

- **P3.4 stands:** ordinary crust <3.88% at M2, <4.33% at O1, tremor exceeds those
  bounds by 5.6× and 3.3×.
- **The band prediction remains untested.** Long-period bounds are 5–8%, far too
  loose to exclude the predicted response.
- **Reaching 1% needs ~400,000 events.** No partition of 25,962 gets there.

---

## 2026-08-22 — P3.2: mechanism-resolved test. Third null on ordinary crust.

**Code:** `eqf::gcmt`, `examples/gcmt_coulomb_sign.rs`. 26 tests.

### The test mechanisms make possible

P3.4 and P3.5 used raw tidal phase, which is **blind to whether the tide loads or
unloads each fault** — a tide promoting failure on a thrust in Japan simultaneously
inhibits it on a normal fault in Greece, and pooled globally those cancel.

GCMT gives 67,263 moment tensors, 18,310 at Mw ≥ 5.5 over 49 years. Resolving ΔCFS
onto each event's *own* fault plane allows a question raw phase cannot pose:

> **Do earthquakes occur preferentially when ΔCFS on their own fault is positive?**

Under no effect, 50%. This aligns the feature with the physics rather than
discarding events, so it is not the stratification P3.5 showed makes bounds worse.

### Result

| | ΔCFS > 0 | null median | p | mean ΔCFS | p |
|---|---|---|---|---|---|
| **Nodal plane 1** | 49.295% | 49.208% | 0.41 | +53.1 Pa | 0.18 |
| **Nodal plane 2** | 48.329% | 48.591% | 0.76 | +42.5 Pa | 0.31 |

**Nothing significant on either plane.**

Worth noting the null medians sit *below* 50% (49.2%, 48.6%). ΔCFS is not
symmetrically distributed across the global population of fault geometries, so the
naive binomial 50% would have been the wrong reference — a fourth-trap-shaped error
the block-shift null caught by construction.

Both nodal planes were run because a moment tensor gives two equally valid planes
and does not say which broke. Neither shows an effect; the ambiguity is not hiding
one.

### Three independent methods, three nulls

| Test | Approach | Result |
|---|---|---|
| P3.4 | Raw tidal phase, global, longitude-corrected | null; <3.88% at M2 |
| P3.5 | Depth-stratified (pre-registered prediction) | null; prediction unsupported |
| **P3.2** | **Mechanism-resolved ΔCFS sign** | **null** |

These are not three slices of one test. They use different features, different
statistics and different nulls. All three agree.

### What the project can now actually say

**Tidal stress measurably modulates slow seismicity, and we cannot detect it in
ordinary earthquakes.**

- **Positive, replicated:** M2, N2 and O1 at Parkfield and Cascadia — two sites
  differing in tectonics, geography, epoch and detection method. Amplitude scales
  faster than linearly (slope 3.56). The M2 solid Earth tide computes to 595 Pa
  against an independently inferred `Aσ₀` of 600 Pa.
- **Bounded null:** ordinary crust below **3.88% at M2** and **4.33% at O1**, where
  tremor shows 21.7% and 14.5% — at least 3–5× weaker, by three methods.
- **Untested:** the band prediction. Long-period bounds of 5–8% cannot exclude the
  predicted response, and reaching 1% needs ~400,000 events against 25,962.

This is consistent with Beeler & Lockner's nucleation argument and with the
literature's genuinely mixed record — Cochran et al.'s factor-3 result came from
shallow thrust faults under large *ocean* tidal loading, which our solid-tide
calculation does not include.

### Stopping position

Per the pre-commitment in doc 16, this is the point to stop rather than keep
searching. Three methods, one direction. **The publishable output is the bound
plus the moonquake and two-site tremor validations**, and that was agreed before
the results were known.

Two things could change it, and neither is a slice:

1. **Ocean tidal loading.** Cochran et al.'s largest effect came from it. Our
   pipeline models solid Earth tide only. Adding it is real physics, not data
   selection — and doc 12 deferred it on a band-prediction argument that the flat
   `R(ω)` measurement has since undermined.
2. **~400,000 events**, which needs regional catalogues at lower Mc with per-region
   completeness control.

---

## 2026-08-22 — Ocean loading: the gap is real and quantified

**Code:** `scripts/setup-spotl.sh`.

SPOTL (Agnew) built and running — Farrell Green's functions, GOT4.7 global ocean
tide model, land-sea database. Four modern-Mac obstacles are documented in the
script: no Fortran compiler, Homebrew gfortran failing to find `libSystem`, K&R C
in `ispand.c` rejected by clang, and `mapcon f` taking **no filenames** (it reads
stdin and writes `lndsea.ind`/`lndsea.bit` itself with `status='new'`).

### Parkfield, M2

```text
ocean loading strain:  3.88 @ 25.2 deg   3.46 @ 145.4 deg   2.70 @ 130.4 deg   (nanostrain)
solid Earth tide strain:                 9.92                                   (nanostrain)
```

**Ocean loading is ~39% of the solid tide in amplitude — and its phases differ
substantially.**

### Why the phase matters more than the amplitude

Amplitude alone would only rescale our bounds. **Phase is what breaks the test.**

Our ordinary-crust null asked whether events cluster in *solid tide* phase. If the
real driver is total tide = solid + ocean, and ocean loading carries a
site-dependent phase offset, then pooling globally against solid-tide phase
**scrambles the signal by construction** — every site contributes at a different
offset and they average away.

That is a precise mechanism by which three null results could all be wrong
together, and it is not something more events would fix.

### And Parkfield is the *mild* case

Parkfield sits ~40 km inland. Ocean loading grows sharply toward the coast, and
Cochran et al.'s factor-3 result came from **shallow thrust faults in subduction
zones** — coastal, where loading routinely exceeds the solid tide. Our global
earthquake catalogue is full of exactly those.

### Status

The tooling works and the magnitude is established. **The remaining work is
running it at scale**, which needs the loading computed per event site — 18,310
locations. `nloadf` convolves a global model per call, so a HEALPix precompute and
interpolation (doc 06 §1, doc 13 §2a) is the route rather than 18,310 invocations.

Until that lands, the honest statement of our result is unchanged but its scope is
now explicit: **we have bounded the response to solid Earth tides, and the
component that produced the literature's largest effect is measured at ~39% of the
solid tide even at an inland site, with a different phase.**

---

## 2026-08-22 — Ocean loading at scale: the null was on an incomplete forcing

Ocean loading strain computed at **all 18,316 GCMT sites** (Mw ≥ 5.5), M2 and O1,
via SPOTL with GOT4.7. A single `nloadf` call is ~25 ms, so direct evaluation beat
a HEALPix precompute — and avoided interpolating a field that varies sharply across
coastlines.

### The numbers

| | median | 90th pct | max |
|---|---|---|---|
| **M2 ocean loading** | **9.23** | **17.90** | **66.3** |
| O1 ocean loading | 2.10 | 5.88 | 24.0 |
| *M2 solid Earth tide, for comparison* | *9.92* | | |

Nanostrain.

**The median ocean loading is 93% of the solid Earth tide. For the top decile of
sites it is 1.8× the solid tide, reaching 6.7× at the extreme.**

This is not a correction term. For a large fraction of the catalogue it is
comparable to or larger than the component we modelled — and it carries a
**site-dependent phase**.

### What this does to P3.4, P3.5 and P3.2

All three tested whether events align with **solid-tide** phase. If the driver is
total tide, and the omitted term is of comparable magnitude with a phase that
varies site to site, then pooling globally against solid-tide phase averages the
signal away **by construction**.

So the correct statement of our result is not what we have been saying:

> ~~Ordinary crust shows no detectable tidal response below ~3.9%.~~
>
> **Ordinary crust shows no detectable response to the *solid Earth tide* below
> ~3.9%. The ocean loading term, comparable in magnitude and differing in phase,
> was not included in that test.**

Three independent methods agreeing does not help here: they shared the same
incomplete forcing. **Agreement between methods is not agreement about physics
when the methods share an input.** That is worth recording as its own lesson —
it is the seventh entry in the trap list in spirit, arriving through a shared
feature rather than a shared statistic.

### Also explains the literature

Cochran et al.'s factor-3 result came from shallow thrust faults in subduction
zones. Those are exactly the sites in our upper decile, where loading is ~2× the
solid tide. A solid-tide-only analysis would find nothing there, and did.

### Remaining work

The data exists; the integration does not. What is needed:

1. Loading strain at event time from the per-site `(amplitude, phase)` pairs and
   the analytic Doodson argument — `ph_core::doodson` already supplies the
   argument, including the longitude correction.
2. Strain → stress under the free-surface condition (plane stress, `σ_zz ≈ 0`),
   giving the horizontal tensor.
3. Add to the solid-tide tensor, resolve on each event's GCMT plane, re-run the
   ΔCFS sign test.

⚠ **One convention must be verified before trusting the result.** SPOTL reports
"phases are local, lags negative". A sign error there flips the ΔCFS sign for every
event and would invert the answer. Validate against SPOTL's own `hartid` time
series at one site, or run both sign conventions as a stated robustness check — the
same treatment the nodal-plane ambiguity already gets.

---

## 2026-08-22 — Loading phase convention pinned before use

The ocean loading data is only usable if SPOTL's reported phase can be mapped onto
our analytic Doodson argument. A sign error there flips ΔCFS for every event; a
smaller offset merely degrades a phase-sensitive test, which is worse because it
would not announce itself.

**Verified against SPOTL's own `hartid` synthesis** rather than assumed
(`examples/verify_loading_phase.rs`).

### Sign

Correlation against a 12-hour `hartid` series: `cos(χ + φ)` gives **+0.86**,
`cos(χ − φ)` gives **+0.32**. Sign convention is **plus**, unambiguously.

### A constant residual, and what it proves

Fitting amplitude and phase over 30-day hourly series at four sites:

| Site | Longitude | Amp fitted | Amp reported | Phase offset |
|---|---|---|---|---|
| Parkfield | −120.2 | 2.7039 | 2.7008 | **10.95°** |
| Japan | +141.0 | 1.3509 | 1.3494 | **10.99°** |
| Europe | +5.0 | 1.4313 | 1.4296 | **10.98°** |
| New Zealand | +175.0 | 21.9386 | 21.9133 | **10.95°** |

**Amplitudes match to 0.1%. The phase offset is constant to 0.04° across 295° of
longitude.**

The constancy is the useful part. **It independently validates the `2λ` longitude
correction** in `doodson::phase_at_longitude` — a longitude-convention error would
make the offset scale with longitude, and it does not. That correction was added
for the global catalogue on theoretical grounds; it is now confirmed empirically at
four sites.

The residual is a definitional difference in the M2 argument reference, not an
error in either implementation. Applied as a calibration constant it raises the
12-hour correlation from 0.86 to **0.93**:

```text
loading(t) = A · cos( χ_local(t) + φ_SPOTL + 10.97° )
```

11° of M2 is 23 minutes — small enough to have been waved through, large enough to
degrade exactly the test we need it for.

### Status

Everything needed for the total-tide test is now in place and validated:

- Ocean loading `(amplitude, phase)` at all 18,316 sites, M2 and O1
- The convention mapping SPOTL's output onto our argument, verified globally
- GCMT mechanisms and both nodal planes
- Solid-tide tensors, elastic calibration, Coulomb projection, block-shift null

**Remaining is mechanical:** reconstruct loading strain at event times, convert
strain to stress under the free-surface condition, add to the solid-tide tensor,
resolve on each plane, re-run the ΔCFS sign test. No unresolved physics or
convention questions stand in the way.
