//! Seismic catalogue ingestion and tidal response measurement.
//!
//! The physics lives upstream in [`ph_core`] — tidal tensors, Coulomb stress,
//! elastic response, analytic constituent phases, and the statistics. This crate
//! holds what only the seismology programme needs: **catalogue ingestion**, and the
//! measurements built on it.
//!
//! # Why these modules are here and not upstream
//!
//! PlanetaryHarmonics is a library serving four downstream projects. Star Seer does
//! not need a USGS earthquake parser; Cosmic Cypher does not need Parkfield LFE
//! families. A module only one consumer needs is application code, not library
//! code — see `docs/14-repo-architecture.md` upstream.
//!
//! # The validation ladder
//!
//! All four catalogues serve one programme, ordered by how strong the effect is:
//!
//! | Rung | Catalogue | Role |
//! |---|---|---|
//! | 1 | [`apollo`] — deep moonquakes | Known answer. Tidal forcing is *dominant*. |
//! | 2 | [`parkfield`], [`cascadia`] — tremor and LFEs | Strong effect, two independent sites. |
//! | 3 | [`comcat`] — global earthquakes | **The actual question.** Ordinary crust. |
//!
//! Rungs 1 and 2 are **controls** for rung 3: tremor has a short `T_a`, so its
//! short-period response is expected and says nothing about ordinary crust.

pub use ph_core;
pub use rustspice_core;

pub mod apollo;
pub mod cascadia;
pub mod comcat;
pub mod parkfield;
