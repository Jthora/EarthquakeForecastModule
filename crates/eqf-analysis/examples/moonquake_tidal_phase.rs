//! A3/A4: assign deep moonquakes a tidal phase from real geometry, then test.
//!
//!     cargo run --release --example moonquake_tidal_phase
//!
//! Requires kernels and the Apollo catalogue — see `scripts/`.
//!
//! The forcing scalar is the largest eigenvalue of the lunar tidal tensor, i.e.
//! peak tidal extension. Because the Moon is tidally locked, Earth sits near the
//! prime meridian and the tensor's *orientation* barely moves; what varies is its
//! *magnitude*, as GM/d³. So maxima of this scalar are lunar perigees, and phase 0
//! means perigee.

use eqf::{apollo};
use ph_core::{field::TidalField, phase::Forcing, stats};
use rustspice_core::{Et, KernelSet};

fn main() -> rustspice_core::Result<()> {
    let mut ks = KernelSet::new();
    for k in [
        "naif0012.tls",
        "de440s.bsp",
        "pck00011.tpc",
        "gm_de440.tpc",
        "moon_pa_de440_200625.bpc",
        "moon_de440_250416.tf",
    ] {
        ks.add_file(format!("kernels/{k}"))?;
    }
    let mut spice = ks.open()?;

    let csv = std::fs::read_to_string("data/apollo/levent.1008weber.csv")
        .expect("run scripts/fetch-apollo.sh first");
    let events = apollo::parse_levent(&csv);
    let dm = apollo::deep_moonquakes(&events, None);
    let days = dm.times();
    println!("deep moonquakes: {}", days.len());

    // Catalogue times are days since 1969-01-01 UTC. Adding 86400 s per day
    // ignores the leap seconds inserted over 1969-1977 (~7 s total). At monthly
    // periods that is 3e-6 of a cycle — far below anything that matters here.
    let epoch0 = spice.parse_time("1969-01-01T00:00:00")?;
    let to_et = |d: f64| Et(epoch0.0 + d * 86400.0);

    // Sample the forcing every 6 hours across the catalogue span, padded so every
    // event falls inside the bracketing maxima.
    let (t0, t1) = (days[0] - 30.0, days[days.len() - 1] + 30.0);
    let step_days = 0.25;
    let n = ((t1 - t0) / step_days) as usize;
    let sample_days: Vec<f64> = (0..n).map(|i| t0 + i as f64 * step_days).collect();
    let epochs: Vec<Et> = sample_days.iter().map(|&d| to_et(d)).collect();

    let lunar = TidalField::on_moon(&mut spice)?;
    let tensors = lunar.tensors(&mut spice, &epochs)?;
    let values: Vec<f64> = tensors.iter().map(|t| t.eigen().0[2]).collect();
    println!("sampled forcing at {} points, {} d step", values.len(), step_days);

    let forcing = Forcing::new(sample_days, values).expect("forcing");
    println!(
        "forcing maxima: {}   mean interval {:.3} d",
        forcing.maxima().len(),
        forcing.mean_period().unwrap()
    );
    println!("  (anomalistic month is 27.555 d — this is a self-check)\n");

    let (phases, dropped) = forcing.phases(&days);
    println!("phased {} events, dropped {dropped}", phases.len());

    println!("\nSchuster test on tidal phase:");
    for s in stats::schuster(&phases, 4) {
        println!(
            "  order {}  D^2/N = {:>8.2}   analytic p = {:.3e}   phase {:>7.1} deg",
            s.order,
            s.d_squared / phases.len() as f64,
            s.p_value,
            s.phase.to_degrees()
        );
    }

    // A4: the time-shift null. Now meaningful, because the forcing is genuinely
    // quasi-periodic — successive perigee intervals differ, so a global shift
    // does not merely rotate the phase cluster.
    let observed = stats::schuster(&phases, 1)[0].d_squared;
    let offsets = stats::shift_offsets(200, 40.0, 1200.0, &[27.555, 205.89], 3.0);
    println!("\ntime-shift null: {} offsets", offsets.len());

    let null = stats::time_shift_null(&days, &offsets, 1, |t| {
        forcing.phase_at(t).unwrap_or(0.0)
    });
    let p = null.p_value(observed);
    println!(
        "  observed D^2 = {:.1}   null max = {:.1}   empirical p = {:.4}",
        observed,
        null.samples.last().copied().unwrap_or(0.0),
        p
    );
    println!(
        "  (floor is 1/(n+1) = {:.4})",
        1.0 / (null.samples.len() as f64 + 1.0)
    );

    Ok(())
}
