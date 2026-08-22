# 17 — Next Steps

Written 2026-08-22, at the stopping point the pre-commitment anticipated. Three
independent methods return null on ordinary crust; tremor is positive and
replicated. This plans what follows.

---

## The one real gap: ocean tidal loading

**Everything we have computed models the solid Earth tide only.**

That is not a detail. **Cochran, Vidale & Tanaka (2004) found their factor-3
effect specifically in shallow thrust faults where *ocean* tidal loading is
large** — the strongest positive result in the literature, driven by the
component we omit. In coastal and subduction settings ocean loading frequently
**exceeds** the solid tide.

So the honest reading of our ordinary-crust null is narrower than it looks:

> We have bounded the response to **solid Earth tides** below 3.88% at M2. We have
> not tested the component that produced the literature's largest effect.

A reviewer would name this immediately, and they would be right to.

### The deferral no longer holds

Doc 12 deferred ocean loading on the argument that the band prediction made
short-period stress unimportant. **The flat `R(ω)` measurement undermined that
argument** — no band limit was found, so there is no longer a reason to think
short-period forcing is irrelevant. The premise for deferring it is gone.

### What it takes

| Step | Notes |
|---|---|
| Ocean tide model | FES2014, TPXO9 or GOT — all free for research |
| Farrell Green's functions | Load Love numbers for an elastic Earth |
| Convolution | Global integral of load × Green's function, per site |
| Validation | Compare against SPOTL or the M3G/Chalmers loading provider at the two tremor sites, where published coefficients exist |

The convolution is the expensive piece and does not scale to 18,310 distinct event
locations naively. **This is where doc 06's spherical-harmonic architecture finally
earns its place:** compute loading coefficients on a HEALPix grid once, then
interpolate. Doc 13 §2a already identifies this as the one surrogate genuinely
worth building.

**Priority: highest.** It is the only remaining step that could change the
ordinary-crust conclusion, and it is missing physics rather than a hope that more
data helps.

---

## The strongest remaining research idea: time-varying sensitivity

`β(x,t)` — tidal sensitivity as a proxy for proximity to failure — has been the
best idea in the project since doc 01. It has never been tested, because until now
there was nowhere with enough signal.

**There is now.** Parkfield and Cascadia both show strong, replicated response.
Beaucé et al. (2023) report tidal sensitivity **rising ~1.5 years before the M7.1
Ridgecrest earthquake**; Ide et al. (2016) report tides modulating the *size*
distribution, not just the rate.

The test: measure `β` in sliding windows at both sites and ask whether it varies
systematically, and whether excursions precede anything.

Two reasons this is worth doing:

1. **It is the precursor claim**, which is the part of this field with real
   forecasting value.
2. **It goes where the signal is.** Our own results say ordinary crust is where the
   effect is undetectable and tremor is where it is 20%.

⚠ Sliding windows shrink the sample per window, and P3.5 established that slicing
degrades bounds. **Pre-register the window length and the excursion criterion
before looking.**

---

## Write-up

The pre-commitment in doc 16 was explicit: the bound plus the validations are
publishable *whatever the result*. That is now due.

**What stands on its own:**

- **Methodological:** six documented traps, each of which produced a confident wrong
  answer — the worst reporting p = 10⁻⁸⁹ where the correct answer was p = 0.70. The
  literature's mixed record on tidal triggering becomes explicable if null choice
  can swing an answer by 88 orders of magnitude. This may be the most useful thing
  the project has produced.
- **Positive:** M2/N2/O1 replicated across two sites differing in tectonics,
  geography, epoch and detection method. Amplitude scaling faster than linear
  (slope 3.56). M2 solid tide computing to 595 Pa against an independently inferred
  `Aσ₀` of 600 Pa.
- **Bounded null:** ordinary crust <3.88% at M2 by three independent methods —
  explicitly for the solid tide only.

**Sequencing:** ocean loading first. It would change the central claim from
"undetectable" to either "undetectable including ocean loading" — much stronger —
or to a positive result. Writing up before it invites the obvious objection.

---

## Deferred, and why: earthquake forecasting

ETAS, the residual model and CSEP evaluation were this repository's original
purpose. They should wait.

**Our own results say the features do not work for ordinary earthquakes.** Fitting
`λ = λ_ETAS · exp(f_θ(features))` on features bounded below 3.88% would be modelling
noise, and the frozen-ETAS design guarantees it would report approximately zero —
correctly, and uninformatively.

**The honest redirect is to apply the forecasting machinery where the signal is.**
Tremor and slow slip show 12–22% modulation, ETS episodes are quasi-periodic, and
`β(t)` is a genuine precursor candidate. Slow-slip forecasting is legitimate science
and we have the measurements to support it.

That is a real change of target and worth deciding deliberately rather than
drifting into. **Flagging it as a decision, not making it.**

---

## Order

```
1. Ocean tidal loading        ← the gap; could change the conclusion
2. Re-run P3.2 and P3.4 with loading included
3. Write up                    ← whichever way step 2 goes
   ‖
   beta(t) at Parkfield and Cascadia   (parallel; independent of 1-3)
```

**If step 2 stays null**, the result is a clean, well-controlled bound on
celestial–terrestrial coupling in ordinary crust, backed by a validated instrument
and a replicated positive in slow seismicity. That was pre-agreed as a publishable
outcome and it remains one.

## Standing rules, unchanged

- **Specify the null before running, and state what structure it preserves.**
- **Stop slicing.** P3.5 established that each stratification loosens the bound.
- Improving a *feature* (ocean loading, mechanisms) is not slicing. Partitioning
  *data* is.
