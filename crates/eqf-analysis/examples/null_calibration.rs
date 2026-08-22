//! Does the block-shift null actually have a 5% false-positive rate?
//!
//!     cargo run --release --example null_calibration
//!
//! # Why this was overdue
//!
//! The block-shift null was invented for this project, argued for on structural
//! grounds, and then relied on for a dozen conclusions — **without ever being
//! calibrated.** Six earlier nulls in this project failed in ways that looked
//! decisive. There is no reason this one is exempt from testing.
//!
//! # Method
//!
//! Generate synthetic catalogues with **no tidal signal whatsoever** but realistic
//! structure, run the full pipeline, and check the p-value distribution.
//!
//! - **Poisson** — no clustering. The easy case.
//! - **Hawkes** — background plus Omori-decaying aftershocks. This is what breaks
//!   naive nulls: strong temporal correlation, non-independent samples.
//! - **Diurnal** — Hawkes plus a 24 h detection modulation. A *non-tidal*
//!   periodicity the null must not mistake for M2.
//!
//! A calibrated null gives **uniform p-values** and rejects at exactly the nominal
//! rate. Anti-conservative (too many false positives) would invalidate results;
//! conservative merely costs power.
//!
//! A power check follows: inject a known modulation and measure detection rate.

use ph_core::doodson;

const N_TRIALS: usize = 200;
const N_EVENTS: usize = 8_000;
const SPAN_DAYS: f64 = 7305.0;
const NULL_SHIFTS: usize = 100;
const LON: f64 = -120.0;

struct Rng(u64);
impl Rng {
    fn new(seed: u64) -> Self {
        Rng(seed | 1)
    }
    fn next_f64(&mut self) -> f64 {
        self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        ((self.0 >> 11) as f64) / ((1u64 << 53) as f64)
    }
}

/// Poisson times, uniform over the span. No structure at all.
fn poisson(n: usize, rng: &mut Rng) -> Vec<f64> {
    let mut t: Vec<f64> = (0..n).map(|_| rng.next_f64() * SPAN_DAYS).collect();
    t.sort_by(|a, b| a.partial_cmp(b).unwrap());
    t
}

/// Hawkes-like: background events each spawning Omori-decaying aftershocks.
/// Produces the heavy temporal clustering that defeats independence-assuming tests.
fn hawkes(n: usize, rng: &mut Rng) -> Vec<f64> {
    let mut t = Vec::with_capacity(n + 64);
    while t.len() < n {
        let parent = rng.next_f64() * SPAN_DAYS;
        t.push(parent);
        // Aftershock count ~ geometric; times ~ Omori p=1.1, c=0.01 d.
        let mut k = 0;
        while rng.next_f64() < 0.6 && k < 40 {
            let u = rng.next_f64().max(1e-12);
            let dt = 0.01 * (u.powf(-1.0 / 0.1) - 1.0);
            if dt < SPAN_DAYS {
                t.push((parent + dt).min(SPAN_DAYS));
            }
            k += 1;
        }
    }
    t.truncate(n);
    t.sort_by(|a, b| a.partial_cmp(b).unwrap());
    t
}

/// Hawkes plus a 24 h detection modulation — a non-tidal periodicity.
fn hawkes_diurnal(n: usize, rng: &mut Rng) -> Vec<f64> {
    let base = hawkes(n * 2, rng);
    // Rejection-sample against a 1 + 0.4 cos(2*pi*t) day-cycle acceptance.
    let mut out = Vec::with_capacity(n);
    for t in base {
        let acc = 0.5 * (1.0 + 0.4 * (std::f64::consts::TAU * t.fract()).cos());
        if rng.next_f64() < acc {
            out.push(t);
        }
        if out.len() == n {
            break;
        }
    }
    out
}

/// Inject a genuine M2 modulation of fractional amplitude `eps`.
fn with_m2(n: usize, eps: f64, rng: &mut Rng) -> Vec<f64> {
    let c = doodson::constituent("M2").unwrap();
    let mut out = Vec::with_capacity(n);
    while out.len() < n {
        let t = rng.next_f64() * SPAN_DAYS;
        let acc = (1.0 + eps * c.phase_at_longitude(t, LON).cos()) / (1.0 + eps);
        if rng.next_f64() < acc {
            out.push(t);
        }
    }
    out.sort_by(|a, b| a.partial_cmp(b).unwrap());
    out
}

fn power_at(c: &doodson::Constituent, days: &[f64], shift: &dyn Fn(f64) -> f64) -> f64 {
    let (mut a, mut b) = (0.0f64, 0.0f64);
    for &t in days {
        let (s, co) = c.phase_at_longitude(t + shift(t), LON).sin_cos();
        a += co;
        b += s;
    }
    (a * a + b * b) / days.len() as f64
}

/// The block-shift p-value, exactly as used throughout the project.
fn block_shift_p(days: &[f64], rng: &mut Rng) -> f64 {
    let c = doodson::constituent("M2").unwrap();
    let period = c.period_days();
    let block = (4.0 * period).max(30.0);
    let (t0, t1) = (days[0], days[days.len() - 1]);
    let nb = (((t1 - t0) / block).floor() as usize) + 2;

    let obs = power_at(c, days, &|_| 0.0);
    let mut ge = 0usize;
    for _ in 0..NULL_SHIFTS {
        let offs: Vec<f64> = (0..nb).map(|_| rng.next_f64() * period).collect();
        let v = power_at(c, days, &|t| {
            offs[(((t - t0) / block).floor() as usize).min(offs.len() - 1)]
        });
        if v >= obs {
            ge += 1;
        }
    }
    (ge as f64 + 1.0) / (NULL_SHIFTS as f64 + 1.0)
}

fn report(label: &str, ps: &[f64]) {
    let n = ps.len() as f64;
    let rate = |a: f64| ps.iter().filter(|&&p| p < a).count() as f64 / n;
    // Kolmogorov-Smirnov against uniform.
    let mut s = ps.to_vec();
    s.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let ks = s
        .iter()
        .enumerate()
        .map(|(i, &p)| {
            let f = (i + 1) as f64 / n;
            (f - p).abs().max((p - i as f64 / n).abs())
        })
        .fold(0.0, f64::max);
    let crit = 1.36 / n.sqrt();
    println!(
        "{label:<22} p<0.01 {:.3}  p<0.05 {:.3}  p<0.10 {:.3}   KS {:.3} (crit {:.3}) {}",
        rate(0.01),
        rate(0.05),
        rate(0.10),
        ks,
        crit,
        if ks < crit { "uniform" } else { "NOT UNIFORM" }
    );
}

fn main() {
    println!(
        "{N_TRIALS} trials, {N_EVENTS} events, {:.0} yr span, {NULL_SHIFTS} shifts\n",
        SPAN_DAYS / 365.25
    );
    println!("=== size: synthetic catalogues with NO tidal signal ===");
    println!("nominal rates are 0.01 / 0.05 / 0.10\n");

    for (label, gen) in [
        ("Poisson", 0u8),
        ("Hawkes (clustered)", 1),
        ("Hawkes + diurnal", 2),
    ] {
        let mut rng = Rng::new(0xCA11B + gen as u64 * 977);
        let ps: Vec<f64> = (0..N_TRIALS)
            .map(|_| {
                let d = match gen {
                    0 => poisson(N_EVENTS, &mut rng),
                    1 => hawkes(N_EVENTS, &mut rng),
                    _ => hawkes_diurnal(N_EVENTS, &mut rng),
                };
                block_shift_p(&d, &mut rng)
            })
            .collect();
        report(label, &ps);
    }

    println!("\n=== power: injected M2 modulation ===");
    for eps in [0.005, 0.01, 0.02, 0.05] {
        let mut rng = Rng::new(0x9051 + (eps * 1e4) as u64);
        let hits = (0..60)
            .filter(|_| block_shift_p(&with_m2(N_EVENTS, eps, &mut rng), &mut rng) < 0.05)
            .count();
        println!("  eps = {:>5.1}%   detected in {hits}/60 trials", eps * 100.0);
    }
}
