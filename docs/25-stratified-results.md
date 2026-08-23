# Stratified analysis — results

Pre-registered in `docs/24-stratified-preregistration.md` (2026-08-23, committed
before the dataset was generated). **Test period never opened.**

---

## 1. Result

All six analyses null, and — unlike the earlier ComCat runs — every one passed
calibration on the first attempt.

| stratum | independent | z sd | max \|z\| | 95th pct | p | Holm |
|---|---|---|---|---|---|---|
| all (reference) | 5,746 | 0.996 | 3.796 | 4.375 | 0.522 | — |
| thrust | 2,195 | 1.043 | 3.827 | 4.409 | 0.502 | 1.00 |
| normal | 2,445 | 0.999 | 3.650 | 4.415 | 0.734 | 1.00 |
| strike-slip | 4,192 | 1.001 | 4.108 | 4.464 | 0.176 | 0.70 |
| shallow (<70 km) | 5,389 | 0.987 | 4.271 | 4.450 | 0.096 | 0.48 |
| deep (≥70 km) | 1,492 | 1.034 | 3.517 | 4.437 | 0.914 | 1.00 |

Smallest raw p is 0.096; Holm-corrected across the five pre-registered strata it
is 0.48. In every stratum the largest of ~9,800 statistics falls below the 95th
percentile of its own permutation null, and in four of six below the *median*.

Gradient-boosted trees on the same strata, against the pre-registered 0.01
bits/event bar:

| stratum | best validation IG |
|---|---|
| thrust | +0.00312 |
| deep | +0.00502 |
| normal | +0.00208 |
| strike-slip | +0.00193 |
| all | +0.00149 |
| shallow | +0.00001 |

Every one below the bar; the largest is half of it, in the stratum declared
underpowered in advance. These sit above the global GBT figure (+0.00006) simply
because smaller validation sets are noisier, which is what the fixed bar is for.

## 2. Bounds

Power was measured, not assumed. Planting a known effect into the thrust stratum
and re-running the whole scan: **β = 0.09 → max\|z\| = 5.546, p = 0.0099,
detected**, against a threshold of 0.089 predicted for that stratum in advance.

| stratum | bound on modulation, per SD |
|---|---|
| all GCMT | 5.5% |
| shallow | 5.7% |
| strike-slip | 6.4% |
| normal | 8.4% |
| thrust | **8.9%** (verified empirically) |
| deep | 10.8% — *declared underpowered in advance* |

So the gap left by the global analysis is now closed down to roughly 6–9% per
mechanism class. An effect had to be large within one class and diluted across
the others to hide from the global test; nothing of that shape is there either.

## 3. Calibration

The whole point of `docs/24` §7. Every stratum's observed z sd landed inside the
permutation's own per-draw spread (quantiles 0.385 to 0.884), and the calendar
check was clean on the GCMT catalogue **even before thinning**:

| thinning | strata | z(date) | z(N–P proxy) | z(day of month) |
|---|---|---|---|---|
| none | 22,722 | −0.866 | +0.677 | +0.566 |
| 500 km / 365 d | 5,746 | −1.115 | −0.109 | +0.296 |

Compare ComCat M4.0+, where the same check gave z(date) = **+3.233** unthinned and
produced a spurious Neptune–Pluto detection at p = 0.005. GCMT's declustered set
is far less spatially clustered, so the design was already well calibrated and
thinning only confirmed it. That is a property of the catalogue, and it is the
reason these results needed no rescuing.

## 4. What this does and does not settle

**Settled.** No chart feature modulates earthquake timing by more than ~6–9% per
standard deviation within any of thrust, normal, strike-slip, or shallow
faulting; nor more than ~4% globally (`docs/23`). Trees find no interaction
reaching a tenth of the pre-registered usefulness bar in any stratum.

**Not settled, and not settleable this way.**

- **Outer-planet configurations.** Still untestable, for the reason in `docs/23`
  §6: annual rates are over-dispersed 2.4× even at M6.0+, so no referent design
  makes a cross-year comparison valid.
- **Effects below ~6%.** The independence ceiling — ~10,000 events globally,
  ~2,000–5,000 per stratum — is a property of seismicity, not of method. More
  catalogue does not buy more independent information.
- **Single-region effects.** A signal confined to one fault system would survive
  all of this. Testing it needs a region whose rate is stationary, which is the
  one avenue left that these designs structurally cannot reach.

## 5. Deviation from pre-registration

Achieved counts are below the census in `docs/24` §2 (e.g. all: 5,746 rather than
7,007) because the census spanned 1976–2024 while the analysis uses only the
pre-test period, 1976–2016. Bounds in §2 above are computed from the achieved
counts, not the projected ones. No stratum was added, dropped, or redefined.
