# Pre-registration: mechanism- and depth-stratified analysis

**Written 2026-08-23, before any stratified test was run.** Committed before the
stratified dataset was generated. Extends `docs/21-preregistration.md`; everything
there that is not overridden here still applies.

---

## 1. Why stratify, and why it is dangerous

`docs/23-ml-results.md` bounds chart features at ~4% per SD **averaged over every
tectonic setting, depth and focal mechanism at once**. An effect confined to one
of them would be diluted by the rest and would not appear.

That is a real gap. It is also exactly how a fishing expedition begins, and how
the six errors in `docs/07-research-log.md` happened. Three commitments:

1. The strata are defined in code (`crates/eqf-dataset/src/strata.rs`) and fixed
   here before any of them is tested.
2. There are **five**, not fifty. Every extra stratum costs a multiplicity penalty
   and, worse, splits the sample.
3. **The power of each stratum is stated below in advance.** A stratum that cannot
   detect what it is looking for contributes nothing but a penalty, so its result
   will be reported as uninformative rather than as a null.

## 2. What this can and cannot find

Detectable effect scales as 1/√n. The global scan needed 10,813 independent
events to reach β = 0.04, which calibrates everything below.

| stratum | declustered | independent | detectable β |
|---|---|---|---|
| all (reference) | 31,011 | 7,007 | 0.050 |
| strike-slip | 15,813 | 5,181 | **0.058** |
| normal | 7,037 | 3,154 | **0.074** |
| thrust | 8,161 | 2,699 | **0.080** |
| shallow (<70 km) | 24,680 | 6,612 | **0.051** |
| deep (≥70 km) | 6,331 | 1,767 | 0.099 |

**The detectable window is narrow and is stated plainly.** These strata see
effects of roughly 6–8% per SD. The global analysis already excludes anything
above 4% on average. So this test is informative only for an effect that is large
*within* a class and diluted *across* classes — an 8% effect in thrust events,
which are 26% of the catalogue, averages to 2% globally and is genuinely
invisible to what has already been run. That is a real gap, but it is a specific
and fairly demanding one, and it is larger than most published tidal-triggering
effects.

**The deep stratum (β = 0.099) is declared underpowered in advance.** Its result
will be reported as uninformative whatever it shows.

## 3. Data

| | |
|---|---|
| catalogue | GCMT `data/gcmt/gcmt.ndk`, 67,263 solutions 1976–2024 |
| magnitude | **no threshold** — see below |
| declustering | Gardner–Knopoff → 31,011 |
| independence thinning | 500 km / 365 days → 7,007 |

No magnitude threshold, deliberately. GCMT is complete globally only at ~M5.5,
and below that the catalogue grows with the network — but the time-stratified
design compares an event only against other dates within its own calendar month,
which conditions out any trend however large. That is what made ComCat M4.0+
usable despite annual counts over-dispersed 294-fold. Here the threshold is a
power choice, not a validity one, and the whole catalogue maximises power.

GCMT is used rather than ComCat because it is the source of the focal mechanisms;
using it for the event list too avoids a cross-catalogue association step and the
matching errors that come with it.

## 4. Strata

Mechanism from the rake of the first nodal plane: within 45° of +90° is thrust,
within 45° of −90° is normal, the rest strike-slip. Rake rather than the Frohlich
P/T/B axis plunges because it is stable under the nodal-plane ambiguity — a thrust
reads near +90° on *both* planes — so no choice between two equally valid planes
is needed. The three classes divide the rake circle into equal quarters/half, so
the definition favours none of them; this is asserted in
`strata::tests::the_three_classes_partition_the_rake_circle_evenly`.

Depth splits at 70 km, the conventional shallow/intermediate boundary.

**Five strata: thrust, normal, strike-slip, shallow, deep.** No cross of mechanism
with depth — that would make nine, and the smallest would be hopeless. No region,
no magnitude band, no plate-boundary class.

## 5. Design, features, metric — unchanged

Time-stratified referents (same calendar month, whole weeks apart). All 9,816
features, no selection. Conditional score test with a permutation null, plus
gradient-boosted trees on the same strata. Information gain in bits per event for
the model; family-wise p from permutation for the scan. Test period stays sealed.

## 6. Success criterion

Per stratum: family-wise p < 0.05 from the permutation null, **Holm-corrected
across the five strata**. Reported alongside that stratum's pre-stated detectable
β, so a null can be read as the bound it is rather than as an absence.

For the tree model: validation IG ≥ 0.01 bits, as before.

## 7. Calibration required before any result is believed

Each stratum must pass the same checks that caught the earlier errors:

- **z(date), z(date²), z(day of month)** consistent with N(0,1)
  (`scripts/calendar_check.py`). The M4.0+ Neptune–Pluto result at p = 0.005 was
  killed by this and nothing else.
- **Observed z sd** within the permutation's per-draw spread. An inflated sd means
  correlated strata, and the result is not reportable until thinning fixes it.
- **Planted-effect recovery** in the stratum, confirming the measured detectable β.

If a stratum fails calibration, no result is reported for it. A miscalibrated
null is how p = 10⁻⁸⁹ happened.

## 8. Declared failure

If every stratum is null, the conclusion recorded is that no chart feature
modulates earthquake timing by more than the stratum-specific bounds in §2, and
that the gap left by the global analysis is now closed down to ~6–8% per class.
That will be written up as plainly as a positive would be.

**A positive in exactly one stratum, with the other four null, is what a false
positive looks like at these multiplicities.** It will be treated as provisional
and reported with the calendar and thinning diagnostics beside it, not as a
finding.
