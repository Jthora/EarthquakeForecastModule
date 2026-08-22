# EarthquakeForecastModule

Probabilistic earthquake forecasting, testing whether celestial, tidal and
gravimetric features improve earthquake-rate estimates over established baselines.

**Start with [HANDOFF.md](HANDOFF.md).** It is written to be read cold and carries
the state of the science, what the upstream library provides, the methodology, and
six documented ways to fool yourself that we have already fallen into.

## Position in the chain

```text
RustSPICE
  └─> PlanetaryHarmonicsModule
        ├─> AstrologyCore ──> Cosmic Cypher, Star Seer, Resonant Finder
        └─> EarthquakeForecastModule          ← this repo
```

The chain **forks** at PlanetaryHarmonics. Forecasting needs tidal tensors and
Coulomb stress; it needs nothing from the interpretive layer and must not depend
on it. The dependency graph is part of the argument.

## Setup

```bash
git submodule update --init --recursive
```

That brings in `PlanetaryHarmonicsModule` and, nested inside it, `RustSPICE`.

Two ways to consume the physics layer — see HANDOFF.md §3:

- **CLI:** `ph-features`, emits CSV with a full provenance header
- **Python:** `maturin build --release` in `modules/PlanetaryHarmonicsModule/crates/ph-py`

## Layout

```text
crates/eqf-analysis/   catalogue ingestion + the measurement programme
  src/                 apollo, parkfield, cascadia, comcat
  examples/            14 analyses, from moonquake validation to the band test
docs/                  research log and methodology
scripts/               data and kernel fetchers -- all public, no credentials
modules/               PlanetaryHarmonicsModule (physics), RustSPICE nested within
```

The physics lives upstream in `ph-core`: tidal tensors, Coulomb stress, elastic
response, analytic constituent phases, statistics. **This repo holds what only
seismology needs** — catalogue parsers and the analyses built on them. A module
only one consumer needs is application code, not library code.

## The validation ladder

| Rung | Catalogue | Role | Result |
|---|---|---|---|
| 1 | Apollo deep moonquakes | known answer | 5/5 periodicities to <0.21% |
| 2 | Parkfield LFEs, Cascadia tremor | strong effect, two sites | M2/N2/O1 significant at both |
| 3 | ComCat global M5.5+ | **ordinary crust — the question** | bounded, inconclusive |

Rungs 1 and 2 are **controls**: tremor has a short `T_a`, so its short-period
response is expected and says nothing about ordinary crust.

## Status

Measurement programme built and run. **No forecasting code yet** — ETAS, the
residual model, β(x,t) and CSEP are still to come (HANDOFF.md §4).

Current position: ordinary crust responds **<3.88% at M2** and **<4.33% at O1**,
where tremor shows 21.7% and 14.5%. The band prediction remains **untested** —
long-period bounds of 5–8% are too loose, and reaching 1% needs roughly 400,000
events against the 25,962 available.

## Reproducing

```bash
git submodule update --init --recursive
./scripts/fetch-kernels.sh && ./scripts/fetch-apollo.sh   # small
./scripts/fetch-parkfield.sh && ./scripts/fetch-cascadia.sh && ./scripts/fetch-comcat.sh
cargo run --release --example moonquake_periodogram
cargo run --release --example band_prediction_test
```

## Terminology

**Forecast, not prediction.** In seismology "prediction" denotes deterministic
time/place/magnitude claims and is largely discredited. This project does
probabilistic rate estimation, which is what CSEP evaluates.
