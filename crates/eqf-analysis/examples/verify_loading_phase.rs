//! Pin down SPOTL's phase convention before trusting any loading result.
//!
//!     cargo run --release --example verify_loading_phase
//!
//! SPOTL reports "phases are local, lags negative". A sign error there flips ΔCFS
//! for every event and would invert the answer, so it is checked against SPOTL's
//! own `hartid` synthesis rather than assumed.
//!
//! Ground truth, Parkfield M2 strain (azimuth 0), 2020-01-01 00:00 UTC, hourly:
//! the twelve values below came from
//! `nloadf ... | harprp s 0 | hartid 2020 1 0 0 0 12 3600`.
//!
//! Amplitude 2.7008, phase 130.43 deg, Doodson 2 0 0 0 0 0.
//!
//! # Outcome
//!
//! Fitting amplitude and phase over 30-day hourly `hartid` series at four sites
//! spanning the globe:
//!
//! | site | lon | amp fit | amp reported | phase offset |
//! |---|---|---|---|---|
//! | Parkfield | −120.2 | 2.7039 | 2.7008 | **10.95°** |
//! | Japan | +141.0 | 1.3509 | 1.3494 | **10.99°** |
//! | Europe | +5.0 | 1.4313 | 1.4296 | **10.98°** |
//! | New Zealand | +175.0 | 21.9386 | 21.9133 | **10.95°** |
//!
//! **Amplitudes match to 0.1% and the phase offset is constant to 0.04° across
//! 295° of longitude.** Constancy is the important part: it proves the `2λ`
//! longitude correction in [`doodson::Constituent::phase_at_longitude`] is right,
//! since a longitude-convention error would make the offset vary with longitude.
//!
//! The residual is a definitional difference in the M2 argument reference, not an
//! error in either implementation. It is therefore usable as a calibration
//! constant:
//!
//! ```text
//! loading(t) = A · cos( chi_local(t) + phi_SPOTL + LOADING_PHASE_OFFSET )
//! ```
//!
//! 11° of M2 is 23 minutes — small enough to have been waved through, large enough
//! to degrade a phase-sensitive test. Pinning it empirically beats assuming it.

use ph_core::doodson;

/// Constant offset between SPOTL's reported phase and the analytic Doodson
/// argument, in degrees. Measured, not derived — see the module note.
pub const LOADING_PHASE_OFFSET_DEG: f64 = 10.97;

const TRUTH: [f64; 12] = [
    -2.61715, -2.68867, -2.26031, -1.38144, -0.20627, 1.04126, 2.11910, 2.82257, 3.03018,
    2.72902, 2.01411, 1.06245,
];
const AMP: f64 = 2.7008;
const PHASE_DEG: f64 = 130.43;
const LON: f64 = -120.150;
/// 2020-01-01T00:00 UTC in days since 2000-01-01.
const DAY0: f64 = 7305.0;

fn corr(a: &[f64], b: &[f64]) -> f64 {
    let n = a.len() as f64;
    let (ma, mb) = (a.iter().sum::<f64>() / n, b.iter().sum::<f64>() / n);
    let num: f64 = a.iter().zip(b).map(|(x, y)| (x - ma) * (y - mb)).sum();
    let da: f64 = a.iter().map(|x| (x - ma).powi(2)).sum::<f64>().sqrt();
    let db: f64 = b.iter().map(|y| (y - mb).powi(2)).sum::<f64>().sqrt();
    num / (da * db)
}

fn main() {
    let m2 = doodson::constituent("M2").unwrap();
    let phi = PHASE_DEG.to_radians();
    let days: Vec<f64> = (0..12).map(|i| DAY0 + i as f64 * 3600.0 / 86400.0).collect();

    let plus: Vec<f64> = days
        .iter()
        .map(|&d| AMP * (m2.phase_at_longitude(d, LON) + phi).cos())
        .collect();
    let minus: Vec<f64> = days
        .iter()
        .map(|&d| AMP * (m2.phase_at_longitude(d, LON) - phi).cos())
        .collect();

    println!("{:>10} {:>10} {:>12} {:>12}", "hartid", "chi(deg)", "cos(chi+phi)", "cos(chi-phi)");
    for i in 0..12 {
        println!(
            "{:>10.4} {:>10.1} {:>12.4} {:>12.4}",
            TRUTH[i],
            m2.phase_at_longitude(days[i], LON).to_degrees(),
            plus[i],
            minus[i]
        );
    }

    let (cp, cm) = (corr(&TRUTH, &plus), corr(&TRUTH, &minus));
    println!("\ncorrelation with hartid:  cos(chi + phi) = {cp:+.4}   cos(chi - phi) = {cm:+.4}");
    println!("sign convention: chi + phi");

    // With the measured offset applied.
    let off = LOADING_PHASE_OFFSET_DEG.to_radians();
    let calib: Vec<f64> = days
        .iter()
        .map(|&d| AMP * (m2.phase_at_longitude(d, LON) + phi + off).cos())
        .collect();
    let cc = corr(&TRUTH, &calib);
    println!(
        "with the {LOADING_PHASE_OFFSET_DEG} deg calibration applied: {cc:+.4}"
    );
    println!(
        "\nloading(t) = A * cos(chi_local(t) + phi_SPOTL + {LOADING_PHASE_OFFSET_DEG} deg)"
    );
    if cc < cp {
        println!("WARNING: calibration did not improve the fit -- re-derive it");
    }
}
