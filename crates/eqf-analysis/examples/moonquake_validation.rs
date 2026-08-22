//! Phase 1 validation — recover known deep moonquake periodicities.
//!
//! Deep moonquakes are tidally driven, and their periods have been known since
//! Lammlein (1977). If this pipeline cannot recover them from the Apollo
//! catalogue, nothing downstream is trustworthy (`docs/12-build-plan.md`).
//!
//! Run with:
//! ```text
//! ./scripts/fetch-apollo.sh
//! cargo run --example moonquake_validation
//! ```

use eqf::apollo;
use ph_core::stats::{self, Power};

/// Periods reported in the deep moonquake literature, in days.
const KNOWN: &[(&str, f64)] = &[
    ("half-month", 13.6),
    ("draconic month", 27.212),
    ("anomalistic month", 27.5546),
    ("synodic month", 29.5306),
    ("solar perturbation", 206.0),
    ("draconic/anomalistic beat", 2190.0),
];

fn main() {
    let path = concat!(env!("CARGO_MANIFEST_DIR"), "/../../data/apollo/levent.1008weber.csv");
    let csv = match std::fs::read_to_string(path) {
        Ok(s) => s,
        Err(e) => {
            eprintln!("cannot read {path}: {e}\nrun ./scripts/fetch-apollo.sh first");
            std::process::exit(1);
        }
    };

    let events = apollo::parse_levent(&csv);
    let all = apollo::deep_moonquakes(&events, None);
    let a1 = apollo::deep_moonquakes(&events, Some(1));

    println!("Apollo PSE catalogue: {} rows parsed", events.len());
    println!("  deep moonquakes (T2 == A): {}", all.len());
    println!("  nest A1:                   {}", a1.len());

    // Derived beats, as a check that the "known" periods are not folklore.
    let (drac, anom, syn) = (27.212, 27.5546, 29.5306);
    let beat = |a: f64, b: f64| 1.0 / (1.0 / a - 1.0 / b).abs();
    println!("\nderived from lunar month lengths:");
    println!("  draconic/anomalistic beat  = {:8.1} d  (literature: ~2190, 6 yr)", beat(drac, anom));
    println!("  anomalistic/synodic beat/2 = {:8.1} d  (literature: ~206)", beat(anom, syn) / 2.0);

    for (label, cat) in [("all nests", &all), ("nest A1", &a1)] {
        let times = cat.times();
        println!("\n=== {label}  (N = {}) ===", times.len());

        let periods = stats::log_periods(5.0, 3000.0, 6000);
        let spectrum = stats::periodogram(&times, &periods, 1);
        let found = stats::peaks(&spectrum, 0.0);

        println!("  strongest peaks:");
        for p in found.iter().take(8) {
            println!("    {:9.3} d   power {:9.2}", p.period, p.power);
        }

        println!("  nearest peak to each known period:");
        for (name, want) in KNOWN {
            match nearest(&found, *want) {
                Some(p) => {
                    let err = (p.period - want) / want * 100.0;
                    let mark = if err.abs() < 2.0 { "OK  " } else { "    " };
                    println!(
                        "    {mark}{name:<26} want {want:8.2} d  got {:8.2} d  ({err:+.2}%)  power {:8.2}",
                        p.period, p.power
                    );
                }
                None => println!("        {name:<26} want {want:8.2} d  no peak"),
            }
        }
    }

    println!(
        "\nNote: peaks are candidates, not significance. Establish that with a\n\
         time-shifted null — and mind that it is degenerate against a single\n\
         exact frequency (see stats.rs)."
    );
}

/// Peak closest to `target` in log-period, within 10%.
fn nearest(peaks: &[Power], target: f64) -> Option<Power> {
    peaks
        .iter()
        .filter(|p| (p.period / target).ln().abs() < 0.10_f64.ln_1p())
        .min_by(|a, b| {
            let (da, db) = ((a.period - target).abs(), (b.period - target).abs());
            da.partial_cmp(&db).unwrap()
        })
        .copied()
}
