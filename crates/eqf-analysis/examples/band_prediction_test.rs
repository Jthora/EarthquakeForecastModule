//! P3.4: the band prediction test — ordinary crust against tremor.
//!
//!     cargo run --release --example band_prediction_test
//!
//! # The claim under test
//!
//! Two timescales bound the responsive band. Nucleation duration `t_n` damps
//! response above `1/t_n`, and Beeler & Lockner extrapolate `t_n ≥ 1 yr` for the
//! San Andreas. Ader's critical period `T_a = 2π Aσ₀/τ̇` is roughly 20–200 yr for
//! ordinary crust. So:
//!
//! ```text
//! predicted responsive band, ORDINARY CRUST:   ~1 year to ~200 years
//! ```
//!
//! Every semidiurnal, diurnal, fortnightly and monthly constituent should be
//! **damped**, with response appearing only at Sa and longer.
//!
//! **Tremor is the control, not the test.** Parkfield and Cascadia have short
//! `T_a`, so their M2/N2/O1 response is expected and says nothing about ordinary
//! crust. The discriminating comparison is *shape*:
//!
//! | Earthquakes vs tremor | Verdict |
//! |---|---|
//! | Same shape (short-period response) | band prediction **refuted** |
//! | Mirror image (long-period only) | band prediction **confirmed** |
//! | Nothing anywhere | underpowered — inconclusive |
//!
//! # Longitude
//!
//! Global events need [`doodson::Constituent::phase_at_longitude`]. Tidal phase is
//! local; for semidiurnal constituents a 180° longitude error is a whole cycle.

use eqf::{cascadia, comcat, parkfield};
use ph_core::{doodson};

const NULL_TRIALS: usize = 400;
const BANDS: &[&str] = &["M2", "N2", "O1", "Q1", "Mf", "Msf", "Mm", "Ssa", "Sa"];

struct Rng(u64);
impl Rng {
    fn next_f64(&mut self) -> f64 {
        self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        ((self.0 >> 11) as f64) / ((1u64 << 53) as f64)
    }
}

/// Observed power, null median, empirical p, and null maximum.
///
/// `lons` is `None` for single-site catalogues (Greenwich phase is fine when every
/// event shares a longitude) and `Some` for the global catalogue.
fn analyse(
    name: &str,
    times: &[f64],
    lons: Option<&[f64]>,
    rng: &mut Rng,
) -> (f64, f64, f64, f64) {
    let c = doodson::constituent(name).unwrap();
    let period = c.period_days();
    let block = (4.0 * period).max(30.0);
    let (t0, t1) = (times[0], times[times.len() - 1]);
    let n_blocks = (((t1 - t0) / block).floor() as usize) + 2;

    let power = |shift: &dyn Fn(f64) -> f64| -> f64 {
        let (mut a, mut b) = (0.0f64, 0.0f64);
        for (i, &t) in times.iter().enumerate() {
            let p = match lons {
                Some(l) => c.phase_at_longitude(t + shift(t), l[i]),
                None => c.phase_at(t + shift(t)),
            };
            let (s, co) = p.sin_cos();
            a += co;
            b += s;
        }
        (a * a + b * b) / times.len() as f64
    };

    let observed = power(&|_| 0.0);
    let mut null: Vec<f64> = (0..NULL_TRIALS)
        .map(|_| {
            let offs: Vec<f64> = (0..n_blocks).map(|_| rng.next_f64() * period).collect();
            power(&|t| offs[(((t - t0) / block).floor() as usize).min(offs.len() - 1)])
        })
        .collect();
    null.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let ge = null.iter().filter(|&&x| x >= observed).count();
    (
        observed,
        null[null.len() / 2],
        (ge as f64 + 1.0) / (NULL_TRIALS as f64 + 1.0),
        null[null.len() - 1],
    )
}

/// Fractional rate modulation implied by a Schuster power.
///
/// For rate `1 + ε cos θ` each event contributes mean `ε/2`, so `D²/N ≈ Nε²/4`.
fn epsilon(d2n: f64, n: usize) -> f64 {
    2.0 * (d2n / n as f64).sqrt()
}

fn main() {
    let pk = parkfield::parse_catalog(
        &std::fs::read_to_string("data/parkfield/LFEcat_Apr2001-Apr2024.csv").expect("parkfield"),
    );
    let cs = cascadia::parse_catalog(
        &std::fs::read_to_string("data/cascadia/cascadia_tremor.csv").expect("cascadia"),
    );
    let eq = comcat::select(
        &comcat::parse_catalog(
            &std::fs::read_to_string("data/comcat/global_m55.csv").expect("comcat"),
        ),
        5.5,
        None,
    );

    let mut pk_t: Vec<f64> = pk.iter().map(|e| e.day).collect();
    pk_t.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let mut cs_t: Vec<f64> = cs.iter().map(|e| e.day).collect();
    cs_t.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let eq_t: Vec<f64> = eq.iter().map(|q| q.day).collect();
    let eq_l: Vec<f64> = eq.iter().map(|q| q.lon_deg).collect();

    let span = |t: &[f64]| (t[t.len() - 1] - t[0]) / 365.25;
    println!("Parkfield LFEs  (tremor, control) : {:>9} events, {:.1} yr", pk_t.len(), span(&pk_t));
    println!("Cascadia tremor (tremor, control) : {:>9} events, {:.1} yr", cs_t.len(), span(&cs_t));
    println!("ComCat M5.5+    (ORDINARY CRUST)  : {:>9} events, {:.1} yr\n", eq_t.len(), span(&eq_t));

    println!(
        "{:<6} {:>8} | {:>8} {:>7} | {:>8} {:>7} | {:>9} {:>7}",
        "band", "period", "PK rat", "PK p", "CS rat", "CS p", "EQ rat", "EQ p"
    );
    println!("{}", "-".repeat(76));

    let mut rng = Rng(0xBA4D);
    let (mut eq_short, mut eq_long) = (0usize, 0usize);
    let mut limits: Vec<(String, f64, f64, f64)> = Vec::new();
    for name in BANDS {
        let (po, pm, pp, _) = analyse(name, &pk_t, None, &mut rng);
        let (co, cm, cp, _) = analyse(name, &cs_t, None, &mut rng);
        let (eo, em, ep, ex) = analyse(name, &eq_t, Some(&eq_l), &mut rng);
        limits.push((name.to_string(), epsilon(ex, eq_t.len()), epsilon(po, pk_t.len()), epsilon(co, cs_t.len())));
        let period = doodson::constituent(name).unwrap().period_days();
        println!(
            "{:<6} {:>8.3} | {:>8.1} {:>6.4}{} | {:>8.1} {:>6.4}{} | {:>9.1} {:>6.4}{}",
            name, period,
            po / pm, pp, if pp < 0.05 { "*" } else { " " },
            co / cm, cp, if cp < 0.05 { "*" } else { " " },
            eo / em, ep, if ep < 0.05 { "*" } else { " " },
        );
        if ep < 0.05 {
            if period < 2.0 { eq_short += 1 } else { eq_long += 1 }
        }
    }

    println!("\nordinary crust: {eq_short} short-period and {eq_long} long-period constituents significant");
    println!(
        "verdict: {}",
        match (eq_short, eq_long) {
            (0, 0) => "nothing significant -- underpowered or no effect; INCONCLUSIVE",
            (0, _) => "long-period only, mirroring the prediction -- CONSISTENT with the band prediction",
            (_, 0) => "short-period only, same shape as tremor -- REFUTES the band prediction",
            _ => "both bands respond -- neither shape; needs interpretation",
        }
    );
    println!("null floor 1/(n+1) = {:.4}", 1.0 / (NULL_TRIALS as f64 + 1.0));

    // A null result is only worth something with a detection limit attached.
    // The smallest modulation we could have called significant is the one that
    // would have exceeded the null maximum.
    println!("\nfractional rate modulation, eps = 2*sqrt((D2/N)/N):");
    println!(
        "{:<6} {:>14} {:>12} {:>12} {:>10}",
        "band", "EQ limit", "PK observed", "CS observed", "PK/EQ"
    );
    println!("{}", "-".repeat(58));
    for (name, eq_lim, pk_eps, cs_eps) in &limits {
        println!(
            "{:<6} {:>13.2}% {:>11.2}% {:>11.2}% {:>10.1}",
            name,
            eq_lim * 100.0,
            pk_eps * 100.0,
            cs_eps * 100.0,
            pk_eps / eq_lim
        );
    }
    println!(
        "\nEQ limit is an UPPER BOUND: ordinary-crust modulation at that constituent\n\
         is below this, or we would have detected it. Where PK/EQ exceeds 1, tremor\n\
         demonstrably responds more strongly than ordinary crust possibly can."
    );
}
