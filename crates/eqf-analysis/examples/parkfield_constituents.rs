//! C1 + D1: LFE response at named tidal constituents, with the artifact gate.
//!
//!     cargo run --release --example parkfield_constituents
//!
//! Rather than a blind peak search, evaluate Schuster power at the **exact**
//! periods of named constituents. That turns doc 08 §13b's validity gate into a
//! direct measurement:
//!
//! - **S2 is exactly 12.000 h**, locked to the day-night cycle. LFE detection is
//!   template matching on continuous data, and its sensitivity varies with
//!   cultural noise, so S2 is the constituent most exposed to a detection
//!   artifact. It also carries the solar *thermal* tide, which is not a body tide.
//! - **M2 is 12.42 h**, so it precesses through local solar time over a lunar
//!   month and decorrelates from any time-of-day artifact.
//! - **S1 (24.000 h) and K1 (23.93 h)** are degenerate with the diurnal artifact
//!   and are effectively unusable. **O1 (25.82 h)** is safe.
//!
//! Power at S1 is therefore a *direct estimate of the artifact floor*: whatever
//! S1 shows is what pure detection bias produces at these sample sizes.

use eqf::{parkfield};
use ph_core::{stats};

/// (name, period in days, usable)
const CONSTITUENTS: &[(&str, f64, bool)] = &[
    ("M2  principal lunar", 0.5175, true),
    ("S2  principal solar", 0.5000, false),
    ("N2  larger elliptic", 0.5274, true),
    ("K1  luni-solar diurnal", 0.9973, false),
    ("S1  solar diurnal", 1.0000, false),
    ("O1  principal lunar diurnal", 1.0758, true),
    ("Mf  lunar fortnightly", 13.661, true),
    ("Msf lunar synodic fortnightly", 14.765, true),
    ("Mm  lunar monthly", 27.555, true),
    ("Ssa solar semiannual", 182.62, true),
    ("Sa  solar annual", 365.26, true),
];

fn main() {
    let path = "data/parkfield/LFEcat_Apr2001-Apr2024.csv";
    let csv = match std::fs::read_to_string(path) {
        Ok(s) => s,
        Err(e) => {
            eprintln!("cannot read {path}: {e}\nrun scripts/fetch-parkfield.sh first");
            std::process::exit(1);
        }
    };

    let events = parkfield::parse_catalog(&csv);
    let fams = parkfield::families(&events);
    let span = {
        let (mut lo, mut hi) = (f64::MAX, f64::MIN);
        for e in &events {
            lo = lo.min(e.day);
            hi = hi.max(e.day);
        }
        hi - lo
    };
    println!(
        "{} events, {} families, span {:.1} d ({:.1} yr)\n",
        events.len(),
        fams.len(),
        span,
        span / 365.25
    );

    // Pooled catalogue first, then the largest family on its own.
    for (label, times) in [
        ("all families pooled", {
            let mut t: Vec<f64> = events.iter().map(|e| e.day).collect();
            t.sort_by(|a, b| a.partial_cmp(b).unwrap());
            t
        }),
        (
            "largest family",
            parkfield::family(&events, &fams[0].0).times(),
        ),
    ] {
        println!("=== {label}: {} events ===", times.len());
        let periods: Vec<f64> = CONSTITUENTS.iter().map(|(_, p, _)| *p).collect();
        let spectrum = stats::periodogram(&times, &periods, 1);

        // The artifact floor, estimated from S1 (exactly 24 h).
        let s1 = spectrum[CONSTITUENTS.iter().position(|c| c.0.starts_with("S1")).unwrap()].power;

        println!(
            "{:<32} {:>9} {:>10} {:>9}",
            "constituent", "period(d)", "power", "vs S1"
        );
        println!("{}", "-".repeat(64));
        for (i, (name, period, usable)) in CONSTITUENTS.iter().enumerate() {
            println!(
                "{:<32} {:>9.4} {:>10.1} {:>9.2}{}",
                name,
                period,
                spectrum[i].power,
                spectrum[i].power / s1,
                if *usable { "" } else { "   [artifact-prone]" }
            );
        }

        let m2 = spectrum[0].power;
        let s2 = spectrum[1].power;
        println!(
            "\nM2/S2 gate: M2 = {m2:.1}, S2 = {s2:.1}, ratio {:.2}",
            m2 / s2
        );
        println!(
            "  {}",
            if m2 > s2 {
                "M2 exceeds S2 -> consistent with a body-tide response"
            } else {
                "S2 exceeds M2 -> suspect detection artifact or thermal tide"
            }
        );
        println!(
            "  artifact floor from S1: {s1:.1}  (power expected ~1 under the null)\n"
        );
    }
}
