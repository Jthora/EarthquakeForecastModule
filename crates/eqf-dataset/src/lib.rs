//! Assembly of the training set: equal-area cells, and matched case-control rows.
//!
//! The feature side lives in `planetary-harmonics-core`; this crate decides *which*
//! (place, time) pairs get featurised and what the negative class means. Those two
//! choices determine whether a result can mean anything, so they are kept separate
//! from the modelling and are tested against a signal-free catalogue.

pub mod cells;
pub mod decluster;
pub mod strata;
pub mod sampling;
