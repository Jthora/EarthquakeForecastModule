# 05 — Research Frontier: Options Outside the Original Scope

Deliberate pass for approaches not in the original project conception. Ordered by
expected value.

---

## 1. Deep moonquakes — the clean-room validation dataset ★

**The best idea in this document.**

Deep moonquakes occur at 700–1200 km depth and are **strongly, uncontroversially
driven by tidal stress** from Earth and the Sun. They cluster at specific
nests and recur with clear monthly (27.55 d), 206-day, and 6-year periodicities.
Unlike terrestrial seismicity, where tidal triggering is a few-percent effect
buried in noise, on the Moon **tidal forcing is the dominant mechanism.**

Data: the **Apollo Passive Seismic Experiment** catalogue (Nakamura et al.),
1969–1977, publicly archived, thousands of classified events.

**Why this matters enormously:**

> If the PlanetaryHarmonics pipeline cannot recover the known deep-moonquake tidal
> periodicities from the Apollo catalogue, the pipeline is broken.

This is a **falsifiable end-to-end test with a known right answer** — something the
terrestrial problem cannot provide. It validates the tidal tensor computation, the
Doodson decomposition, the harmonic encoding, and the point-process machinery all
at once, on data where the effect is large enough to see clearly.

It also has a second benefit: it is a completely defensible publication in its own
right, establishing methodological credibility *before* making any terrestrial
claim. Recommend doing this **first**, before the USGS work.

No zodiac, no interpretation, no controversy — just a validated instrument.

---

## 2. Tectonic tremor and slow slip — high-SNR terrestrial testbed ★

Nonvolcanic tremor is **dramatically more tidally sensitive than ordinary
earthquakes**.

- Rubinstein et al. (2008), *Science* 319 — *Tidal Modulation of Nonvolcanic
  Tremor*: clear tremor pulsing at **12.4 h and 24–25 h**, matching principal
  lunar and lunisolar tides. Small tidal stresses "influence the genesis of tremor
  much more effectively than they do the genesis of normal earthquakes."
- Thomas et al. (2013), *GRL* 40 — tidal triggering of rapid tremor reversals in
  northern Cascadia
- Hawthorne & Rubin (2010), *JGR* 115 — tidal modulation of slow slip in Cascadia
- Ide et al. note tremor rate increases **exponentially** with tidal stress — a
  strongly nonlinear response

**Why it matters:** tremor gives a terrestrial dataset where the effect is large
enough for rapid iteration. Develop and tune the feature pipeline here, then
transfer to the low-SNR earthquake problem. Catalogues exist for Cascadia,
Nankai, and Parkfield (LFE catalogues).

Combined with the moonquake test, this gives a **graded difficulty ladder**:
moonquakes (effect dominant) → tremor (effect strong) → earthquakes (effect small).
That is exactly the right development order and a compelling narrative for a paper.

---

## 3. Induced seismicity as a controlled experiment

Injection-induced seismicity (Oklahoma, Groningen, geothermal fields) offers
something rare: **known, recorded stressing history**. Injection volumes and
pressures are documented.

This gives a semi-controlled test of whether the pipeline correctly separates
a known anthropogenic driver from celestial covariates. If the model attributes
injection-driven seismicity to planetary features, that is a decisive
methodological failure caught cheaply.

Also: these catalogues are dense, well-located, shallow, and have high pore
pressure — conditions under which tidal sensitivity is expected to be elevated.

---

## 4. Hydrological and seasonal loading

Real, published, and usually omitted from celestial-correlation studies — which
makes it a serious **confounder**:

- Seasonal groundwater and snow loading modulate California seismicity
  (Johnson et al.)
- Atmospheric pressure loading
- Reservoir-induced seismicity

**Critical point:** these are *annual* signals. So are **Sa** and **Ssa** tidal
constituents. Any annual-period "celestial" correlation is confounded with
hydrology by default. Hydrological loading must be an explicit covariate, not
ignored — otherwise the strongest apparent celestial signal in the data may simply
be groundwater.

This is probably the most dangerous unaddressed confounder in the project.

---

## 5. Volcanic and geothermal systems

Volcanic seismicity shows elevated tidal sensitivity (magma chambers are
compliant, pore pressure is high). Scholz et al. 2019 explains mid-ocean ridge
triggering through magma chamber inflation.

Useful as another moderate-SNR testbed, and relevant to eruption forecasting as a
downstream application beyond earthquakes.

---

## 6. Mars — InSight marsquakes

The InSight seismometer catalogue (2019–2022) provides a second planetary body
with tidal forcing from Phobos, Deimos, and the Sun, no oceans, and no
hydrological cycle. Some claimed tidal correlations; contested and low event
counts.

Lower priority than the Moon — fewer events, weaker forcing — but a natural third
validation body, and the absence of oceans and hydrology removes two confounders
entirely.

---

## 7. Earth rotation, polar motion, pole tide

- **Pole tide** from Chandler wobble (~433 d) — real centrifugal-potential stress
- **LOD variations** — Bendick & Bilham (2017) claim a 32-year cycle in M7+ rate
  with a ~5-year lead
- IERS provides all series

Physically real (Tier A) with contested interpretation (Tier C). Cheap to include.
The ~433 d Chandler period is close enough to annual and semiannual terms that
careful separation is needed.

---

## 8. b-value / magnitude-distribution modulation

Following Ide et al. (2016): tides may modulate the **size distribution** rather
than the rate. Most studies test rate only.

This is under-explored, requires no new data, and doubles the hypothesis space in
a physically motivated direction. Implement as a magnitude distribution
conditioned on tidal amplitude within the point-process likelihood.

---

## 9. Approaches to treat with caution

**Solar dynamo / planetary tidal forcing of solar activity.** Live literature
(Scafetta; Abreu et al. 2012, *A&A*), with published rebuttals ("Tidally
synchronized solar dynamo: a rebuttal", arXiv:2206.14809). If used at all, cite
both sides and mark Tier C. The proposed Earth-coupling chain is long and weak.

**Geomagnetic / space weather coupling to seismicity.** Claimed correlations exist
but mechanism is poorly established. Tier C at best.

Neither should appear in a first publication.

---

## Recommended sequencing

```
Phase 1  Deep moonquakes  — validate the instrument on a known answer
Phase 2  Tremor / slow slip — tune features on high-SNR terrestrial data
Phase 3  USGS earthquakes  — the hard problem, with a validated pipeline
```

Each phase produces a defensible result. Each de-risks the next. And by Phase 3
the methodology has been validated twice on data where the answer is already
known — which is the strongest possible position from which to make a
terrestrial claim.
