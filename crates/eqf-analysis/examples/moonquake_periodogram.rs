//! Phase 1 validation: recover known deep moonquake periodicities.
//!
//! The Apollo deep moonquake catalogue is the project's known-answer test. Tidal
//! forcing there is *dominant* rather than marginal, and the periodicities have
//! been published since Lammlein (1977). If this pipeline cannot recover them,
//! nothing downstream is trustworthy.
//!
//!     cargo run --example moonquake_periodogram
//!
//! Requires `data/apollo/levent.1008weber.csv` — see `scripts/fetch-apollo.sh`.

use eqf::{apollo};
use ph_core::{stats};

/// Periodicities reported for deep moonquakes, in days.
const KNOWN: &[(&str, f64)] = &[
    ("half-month", 13.606),
    ("draconic month", 27.212),
    ("anomalistic month", 27.555),
    ("synodic month", 29.531),
    ("half full-moon cycle", 205.89),
];

fn main() {
    let path = "data/apollo/levent.1008weber.csv";
    let csv = match std::fs::read_to_string(path) {
        Ok(s) => s,
        Err(e) => {
            eprintln!("cannot read {path}: {e}\nrun scripts/fetch-apollo.sh first");
            std::process::exit(1);
        }
    };

    let events = apollo::parse_levent(&csv);
    let dm = apollo::deep_moonquakes(&events, None);
    let times = dm.times();
    let span = times.last().unwrap() - times.first().unwrap();

    println!("parsed {} catalogue rows", events.len());
    println!("deep moonquakes (T2 == A): {}", dm.len());
    println!("span: {:.1} days ({:.1} years)\n", span, span / 365.25);

    let periods = stats::log_periods(5.0, 3000.0, 20_000);
    let spectrum = stats::periodogram(&times, &periods, 1);
    let peaks = stats::peaks(&spectrum, 20.0);

    println!("{:<12} {:>10} {}", "period (d)", "power", "match");
    println!("{}", "-".repeat(46));
    for p in peaks.iter().take(12) {
        let hit = KNOWN
            .iter()
            .find(|(_, k)| (p.period - k).abs() / k < 0.01)
            .map(|(name, k)| format!("{name} ({k} d)"))
            .unwrap_or_default();
        println!("{:<12.3} {:>10.1} {}", p.period, p.power, hit);
    }

    println!("\nrecovery of known periods:");
    let mut recovered = 0;
    for (name, known) in KNOWN {
        match peaks
            .iter()
            .filter(|p| (p.period - known).abs() / known < 0.01)
            .max_by(|a, b| a.power.partial_cmp(&b.power).unwrap())
        {
            Some(p) => {
                let err = 100.0 * (p.period - known) / known;
                println!(
                    "  {:<24} {:>9.3} d  ({:+.2}%)  power {:.0}",
                    name, p.period, err, p.power
                );
                recovered += 1;
            }
            None => println!("  {name:<24}   NOT RECOVERED"),
        }
    }
    println!("\n{recovered}/{} known periods recovered", KNOWN.len());

    // Long-period peaks are not necessarily physical. The catalogue spans only
    // ~8 years with large observational gaps, and those produce spectral
    // structure of their own. Fewer than ~4 cycles in the record is not enough
    // to distinguish a periodicity from a trend or a gap pattern.
    let suspects: Vec<_> = peaks
        .iter()
        .take(12)
        .filter(|p| p.period > span / 4.0)
        .collect();
    if !suspects.is_empty() {
        println!(
            "\nsuspect long-period peaks (< 4 cycles in a {:.0} d span — \
             likely catalogue artifacts, resolve with the time-shift null):",
            span
        );
        for p in suspects {
            println!(
                "  {:>9.1} d   {:.1} cycles   power {:.0}",
                p.period,
                span / p.period,
                p.power
            );
        }
    }
}
