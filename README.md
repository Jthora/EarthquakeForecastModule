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

## Status

Scaffolding. No forecasting code yet. The immediate work is HANDOFF.md §4.

## Terminology

**Forecast, not prediction.** In seismology "prediction" denotes deterministic
time/place/magnitude claims and is largely discredited. This project does
probabilistic rate estimation, which is what CSEP evaluates.
