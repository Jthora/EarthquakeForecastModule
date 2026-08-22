//! Is the tremor tidal signal real triggering, or detection-threshold bias?
//!
//!     cargo run --release --example detection_bias_test
//!
//! # The alternative explanation
//!
//! Ocean loading deforms the ground and modulates microseism noise. **Custodio et
//! al. (2003, GRL) measured tidal modulation of seismic noise with "the main
//! periodicity coinciding with the tidal component M2"** — the exact frequency
//! where we report a result. Detection threshold tracks noise, so detection
//! capability itself oscillates at M2, producing apparent rate modulation with
//! **zero triggering**.
//!
//! This explains our pattern uncomfortably well:
//!
//! | Catalogue | Near detection threshold? | Our M2 result |
//! |---|---|---|
//! | Parkfield LFEs | yes — template matching at threshold | 21.7% |
//! | Cascadia tremor | yes — envelope cross-correlation | 12.7% |
//! | Global M5.5+ | no — complete | nothing |
//!
//! It also undermines C5's "independent replication": both sites are
//! threshold-limited, so a shared artifact mechanism survives every difference in
//! tectonics, geography and epoch that made the replication look convincing.
//!
//! # Pre-registered prediction
//!
//! Stratify by detection strength — Parkfield by `ccsum` (summed cross-correlation
//! across channels), Cascadia by `num_stas` and magnitude.
//!
//! | Outcome | Reading |
//! |---|---|
//! | Modulation **falls** as detection strength rises | **detection bias** |
//! | Modulation roughly **flat** across strata | **real triggering** |
//! | Modulation **rises** with strength | real, and consistent with Ide's b-value result |
//!
//! Fixed before running. The strongest-detected events are far above threshold, so
//! a detection artifact must vanish there while triggering need not.

use eqf::{cascadia, parkfield};
use ph_core::doodson;

const STRATA: usize = 5;
const NULL_TRIALS: usize = 200;

struct Rng(u64);
impl Rng {
    fn next_f64(&mut self) -> f64 {
        self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        ((self.0 >> 11) as f64) / ((1u64 << 53) as f64)
    }
}

fn power(c: &doodson::Constituent, days: &[f64], lon: f64, shift: &dyn Fn(f64) -> f64) -> f64 {
    let (mut a, mut b) = (0.0f64, 0.0f64);
    for &t in days {
        let (s, co) = c.phase_at_longitude(t + shift(t), lon).sin_cos();
        a += co;
        b += s;
    }
    (a * a + b * b) / days.len() as f64
}

/// Fractional rate modulation and its empirical p, for one stratum.
fn analyse(name: &str, days: &[f64], lon: f64, rng: &mut Rng) -> (f64, f64) {
    let c = doodson::constituent(name).unwrap();
    let period = c.period_days();
    let block = (4.0 * period).max(30.0);
    let (t0, t1) = (days[0], days[days.len() - 1]);
    let nb = (((t1 - t0) / block).floor() as usize) + 2;

    let obs = power(c, days, lon, &|_| 0.0);
    let mut null: Vec<f64> = (0..NULL_TRIALS)
        .map(|_| {
            let offs: Vec<f64> = (0..nb).map(|_| rng.next_f64() * period).collect();
            power(c, days, lon, &|t| {
                offs[(((t - t0) / block).floor() as usize).min(offs.len() - 1)]
            })
        })
        .collect();
    null.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let ge = null.iter().filter(|&&x| x >= obs).count();
    (
        2.0 * (obs / days.len() as f64).sqrt(),
        (ge as f64 + 1.0) / (NULL_TRIALS as f64 + 1.0),
    )
}

/// Sort by `key`, split into equal-count strata, report modulation in each.
fn stratify(label: &str, mut rows: Vec<(f64, f64)>, lon: f64, key: &str, rng: &mut Rng) {
    rows.retain(|(k, _)| k.is_finite());
    if rows.len() < STRATA * 500 {
        println!("{label}: too few usable rows ({})", rows.len());
        return;
    }
    rows.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    let per = rows.len() / STRATA;

    println!("\n=== {label}, stratified by {key} (weakest detection first) ===");
    println!("{:>7} {:>9} {:>12} {:>10} {:>8}", "stratum", "events", key, "M2 eps", "p");
    let mut eps = Vec::new();
    for s in 0..STRATA {
        let hi = if s + 1 == STRATA { rows.len() } else { (s + 1) * per };
        let slice = &rows[s * per..hi];
        let mean_k = slice.iter().map(|(k, _)| k).sum::<f64>() / slice.len() as f64;
        let mut days: Vec<f64> = slice.iter().map(|(_, d)| *d).collect();
        days.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let (e, p) = analyse("M2", &days, lon, rng);
        println!(
            "{:>7} {:>9} {:>12.3} {:>9.2}% {:>8.4}{}",
            s + 1,
            slice.len(),
            mean_k,
            e * 100.0,
            p,
            if p < 0.05 { " *" } else { "" }
        );
        eps.push(e);
    }
    let (first, last) = (eps[0], eps[STRATA - 1]);
    println!(
        "  weakest {:.2}%  ->  strongest {:.2}%   ratio {:.2}",
        first * 100.0,
        last * 100.0,
        last / first
    );
    println!(
        "  {}",
        if last / first < 0.5 {
            "FALLS with detection strength -> consistent with DETECTION BIAS"
        } else if last / first > 1.5 {
            "RISES with detection strength -> real, and not a threshold artifact"
        } else {
            "roughly FLAT across strata -> consistent with real triggering"
        }
    );
}

fn main() {
    let mut rng = Rng(0xDE7EC7);

    let pk = parkfield::parse_catalog(
        &std::fs::read_to_string("data/parkfield/LFEcat_Apr2001-Apr2024.csv").expect("parkfield"),
    );
    stratify(
        "Parkfield LFEs",
        pk.iter().map(|e| (e.ccsum, e.day)).collect(),
        -120.150,
        "ccsum",
        &mut rng,
    );

    let cs = cascadia::parse_catalog(
        &std::fs::read_to_string("data/cascadia/cascadia_tremor.csv").expect("cascadia"),
    );
    stratify(
        "Cascadia tremor",
        cs.iter().map(|e| (e.num_stas, e.day)).collect(),
        -123.0,
        "num_stas",
        &mut rng,
    );
    stratify(
        "Cascadia tremor",
        cs.iter().map(|e| (e.magnitude, e.day)).collect(),
        -123.0,
        "magnitude",
        &mut rng,
    );
}
