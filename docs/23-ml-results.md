# Chart features as earthquake-timing predictors — results

Pre-registered in `docs/21-preregistration.md` (2026-08-22, committed before any
fit). Model grid in `docs/22-model-grid.md`. **The sealed 2017–2024 test period
was never opened** — see §8 for why.

---

## 1. Executive summary

Nine thousand eight hundred and sixteen astrological features, computed with no
selection or physical filtering, tested against 49 years of global seismicity.

**Sub-monthly configurations** — Moon, Sun, Mercury, hour angles, everything that
moves within a calendar month — are testable, and nothing was found. The bound is
quantitative: **no feature modulates earthquake timing by more than about 4% per
standard deviation.** That is comparable to the 3.88% M2 tidal bound this
programme reached by independent means (`docs/19-results.md`).

**Slow configurations** — outer-planet aspects, the substance of classical
astrology — turn out to be **untestable against earthquake catalogues by any
referent design**, for a reason that is about the Earth and not about the method
(§6). That is a weaker statement than "no effect", and the difference matters.

Three separate bugs and one design flaw were found along the way, each of which
had been silently producing a plausible null. Two of them would have been
reported as results.

## 2. What was built

9,816 features per (epoch, site) from `planetary-harmonics-core`: aspects (7,920
= 55 pairs × 24 harmonics × cos/sin × 3 frames), declination parallels (495),
galactic-centre aspects (240), CosmicCypher resonance (216), motion and stations
(245), lunar synodic/anomalistic/draconic phases and node (86), chart-shape
statistics (27), eclipse geometry (7), site-local sidereal time / ascendant /
midheaven / hour angles / altitudes (580). Geocentric, heliocentric, barycentric.

Every feature is a cos/sin pair or a rotation-invariant statistic, so no result
can depend on where the zodiac's zero point is placed.

Validated against external standards rather than internal consistency: GMST
reproduces the defined 18h41m50.548s at J2000; the Greenwich-noon Sun-to-MC
offset is 0.771°, which is the equation of time for 1 January; the lunar node
regresses at −19.34°/yr; 893 planetary stations are found in 49 years where
orbital periods predict ~891.

## 3. The effective sample size is ~10,000, not 488,000

| | |
|---|---|
| ComCat 1976–2024 | 488,215 events |
| after Gardner–Knopoff declustering | 154,302 |
| after thinning to 500 km / 365 days | **10,813** |

Windowed declustering removes aftershock sequences but leaves long-range regional
correlation. Two events months apart in one region share nearly identical values
for any slowly-moving feature, and a within-stratum permutation null cannot see
that because it treats strata as independent.

Going from M5.5+ to M4.0+ multiplies the catalogue twelvefold and adds almost no
independent information. **That is the statistical ceiling on this question, and
it is a property of seismicity.**

## 4. What the calibration diagnostics caught

Each of these produced a believable null before it was found.

**Conditioning (Model A).** Standardisation folded into the coefficients —
algebraically correct, since the mean cancels inside the stratum softmax. It left
L-BFGS facing columns whose variances span 0.5 to 4×10⁸, and at λ=1e-4, where
9,781 features against 8,414 strata should overfit to nearly log₂5 bits, it
reported convergence at 0.0005. *The synthetic null missed it because those
features were all N(0,1) — a calibration set better-conditioned than the real
data catches logical faults, not numerical ones.*

**Float32 centring.** An outer-planet aspect moves less than float32 epsilon
across a 5-day stratum, so its centred residual is quantisation noise and z
explodes. Permutation median max|z| was 21.5; in float64, 3.8.

**Non-exchangeable referents.** With controls at *t* ± k days the controls
straddle the case, so the case sits at its own stratum's temporal centroid and is
closer to the stratum mean than a random member. Observed z sd came out at 0.628
where the permutation null said 1.001 ± 0.032 — twelve standard deviations narrow.
Conservative, so the earlier nulls survived, but the null was wrong. Replaced by
time-stratified referents (same calendar month, whole weeks apart), where
exchangeability holds by construction.

**Clipped referent windows.** Year blocks anchored at 1900 meant cases in 1976–77
kept only referents later than themselves. One-sided sets, and z(date) = −6.4.

## 5. The positive result that wasn't

The corrected M4.0+ design produced a hit: max|z| = 5.525, family-wise p = 0.0050,
and the top thirty features were all Neptune–Pluto aspects, harmonics 5 to 18.

Neptune–Pluto is the slowest pair in the solar system: 0.006°/day. Any feature
built from it is essentially a smooth function of the date. Computing the same
conditional score for plain functions of date, using row metadata only:

| statistic | z | p |
|---|---|---|
| calendar date | +3.233 | 0.0020 |
| date² | −3.719 | 0.0015 |
| **day of month** | **+3.262** | **0.0015** |
| Neptune–Pluto proxy | +3.709 | 0.0015 |

Day of the month scored as high as the astronomy. Under thinning, z(date) fell
from +3.233 to +0.042 and the Neptune–Pluto proxy from +3.709 to −0.778. The
whole effect was correlated strata; variance inflation was 1.43×, and it lands
hardest on the slowest features, which is why the slowest pair won.

## 6. Why slow configurations cannot be tested

Testing an outer-planet aspect means comparing across years, because that is the
only timescale on which it moves. That is valid only if the event was equally
likely to have fallen in any of those years. Declustered annual counts say
otherwise:

| threshold | mean/yr | variance/mean | r(year) |
|---|---|---|---|
| M4.0+ | 3,149 | **294×** | +0.942 |
| M5.0+ | 693 | 8.9× | +0.026 |
| M5.5+ | 248 | 2.6× | +0.273 |
| M6.0+ | 89 | 2.4× | +0.352 |

Even at M6.0+ the annual rate is 2.4× over-dispersed: great earthquakes drive
regional sequences lasting years. No referent design controls a rate that varies
across the very years being compared, and thinning cannot help because the
variation is global — it moves every cell at once. (The M4.0+ figure of 294× with
r = 0.94 is network completeness, which is why the within-month design is correct
at that threshold.)

**Outer-planet configurations are therefore beyond the reach of this catalogue.**

## 7. The calibrated results

### Score scan — full design power, one feature at a time

| dataset | strata | max \|z\| | null median | 95th pct | p |
|---|---|---|---|---|---|
| M5.5+ day-offset | 10,154 | 3.514 | 3.791 | 4.398 | 0.90 |
| M4.0+ day-offset | 119,620 | 3.474 | 3.885 | 4.571 | 0.92 |
| **M4.0+ time-stratified, thinned** | **10,813** | **3.653** | **3.831** | **4.376** | **0.7508** |

In every calibrated run the largest of ~9,800 statistics sits *below the median*
of its own null.

### Power — measured, not assumed

Planting a known effect into the real feature matrix and re-running the whole
scan:

| planted β | max \|z\| | p | |
|---|---|---|---|
| 0.02 | 3.474 | 0.97 | not detected |
| 0.03 | 3.455 | 0.98 | not detected |
| **0.04** | **4.734** | **0.020** | **detected** |
| 0.05 | 5.721 | 0.010 | detected |
| 0.20 | 20.966 | 0.010 | detected |

Detection threshold ≈ **0.035–0.04 log-odds per SD**. Observed: 3.653. So the
result is a bound, not an absence of evidence.

### Model B — gradient-boosted trees, for interactions

| | best validation IG |
|---|---|
| M4.0+ time-stratified, thinned, 8,736 train strata | **+0.00006 bits/event** |

Pre-registered meaningfulness bar: 0.01 bits. This is 170× below it, with every
configuration early-stopping after 1–6 trees. Its capability is separately
established: with a planted β=0.20 and only 1,500 strata it reaches +0.0128 bits,
and with β=0.80, +0.313.

### Model A — conditional logistic

Negative at every λ, and **uninformative**: the same pipeline cannot recover a
planted β=0.50. With 9,781 features against 8,414 strata, L2 must shrink so hard
that real signal goes with the noise. Reported for completeness only.

## 8. The sealed test period was not opened

Pre-registration §8 requires test IG ≥ 0.01 bits with p < 0.01. The best
validation figure across every model and design is +0.00006 bits. The criterion
is unreachable, and scoring the test period would add no information while
spending a resource that can only be spent once.

This is a deliberate abstention and is recorded as a deviation. The 2017–2024
data remains available for a future design — a better one, or a larger catalogue,
or a different question.

## 9. What would change this

- **A catalogue with more independent events.** The binding constraint is
  ~10,000, not 488,000. Only time, or a fundamentally better way to establish
  independence, moves it.
- **Regional analysis where the local rate is stationary.** The global test cannot
  see an effect confined to one tectonic setting.
- **Slow configurations remain open**, not refuted. Testing them needs a hazard
  model that accounts for multi-year rate variation, which is a research problem
  in its own right.

## 10. Multiple comparisons

All pre-registered analyses are null, so Holm correction changes nothing. Stated
for completeness: 3 designs × 3 magnitude thresholds × 2 model classes, with the
headline being M4.0+ time-stratified thinned, score scan and GBT.
