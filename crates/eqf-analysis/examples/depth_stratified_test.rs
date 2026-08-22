//! P3.5: does ordinary-crust response strengthen for shallow events?
//!
//!     cargo run --release --example depth_stratified_test
//!
//! # A single pre-registered prediction, not a scan
//!
//! P3.4 bounded the global M5.5+ response at <3.9% (M2) with nothing significant.
//! That catalogue mixes every depth and fault geometry, which dilutes any effect
//! concentrated in one setting.
//!
//! **Métivier et al. (2009)**, on 442,412 NEIC events, report the tidal phase
//! anomaly is **larger for shallower earthquakes**. **Cochran et al. (2004)** find
//! their factor-3 effect specifically in *shallow thrust* faults where ocean tidal
//! loading is large.
//!
//! So the prediction, fixed before running: **shallow events respond more strongly
//! than deep ones.** One comparison, one split at 70 km (the conventional
//! shallow/intermediate boundary), no scanning over cut depths.
//!
//! If shallow response exceeds deep, that is a positive result on ordinary crust.
//! If not, the bound tightens for the subset where the effect is *predicted* to be
//! strongest, which is more informative than the mixed bound.
//!
//! ⚠ Splitting halves the sample, so each bound is looser by √2. A null here is
//! weaker evidence than P3.4's, not stronger.

use eqf::{comcat};
use ph_core::{doodson};

const NULL_TRIALS: usize = 400;
const BANDS: &[&str] = &["M2", "N2", "O1", "Mf", "Mm", "Sa"];
const SHALLOW_KM: f64 = 70.0;

struct Rng(u64);
impl Rng {
    fn next_f64(&mut self) -> f64 {
        self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        ((self.0 >> 11) as f64) / ((1u64 << 53) as f64)
    }
}

/// Observed power, empirical p, and the null maximum.
fn analyse(name: &str, days: &[f64], lons: &[f64], rng: &mut Rng) -> (f64, f64, f64) {
    let c = doodson::constituent(name).unwrap();
    let period = c.period_days();
    let block = (4.0 * period).max(30.0);
    let (t0, t1) = (days[0], days[days.len() - 1]);
    let n_blocks = (((t1 - t0) / block).floor() as usize) + 2;

    let power = |shift: &dyn Fn(f64) -> f64| -> f64 {
        let (mut a, mut b) = (0.0f64, 0.0f64);
        for (i, &t) in days.iter().enumerate() {
            let (s, co) = c.phase_at_longitude(t + shift(t), lons[i]).sin_cos();
            a += co;
            b += s;
        }
        (a * a + b * b) / days.len() as f64
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
        (ge as f64 + 1.0) / (NULL_TRIALS as f64 + 1.0),
        null[null.len() - 1],
    )
}

fn eps(d2n: f64, n: usize) -> f64 {
    2.0 * (d2n / n as f64).sqrt()
}

fn main() {
    let all = comcat::select(
        &comcat::parse_catalog(
            &std::fs::read_to_string("data/comcat/global_m55.csv").expect("comcat"),
        ),
        5.5,
        None,
    );

    let shallow: Vec<_> = all.iter().filter(|q| q.depth_km <= SHALLOW_KM).collect();
    let deep: Vec<_> = all.iter().filter(|q| q.depth_km > SHALLOW_KM).collect();
    println!(
        "global M5.5+: {} events -> {} shallow (<= {SHALLOW_KM} km), {} deep",
        all.len(),
        shallow.len(),
        deep.len()
    );
    println!("prediction (Metivier 2009): shallow responds more strongly\n");

    let split = |v: &[&comcat::Quake]| -> (Vec<f64>, Vec<f64>) {
        (v.iter().map(|q| q.day).collect(), v.iter().map(|q| q.lon_deg).collect())
    };
    let (sd, sl) = split(&shallow);
    let (dd, dl) = split(&deep);

    println!(
        "{:<6} {:>8} | {:>9} {:>8} {:>7} | {:>9} {:>8} {:>7}",
        "band", "period", "shal eps", "shal lim", "p", "deep eps", "deep lim", "p"
    );
    println!("{}", "-".repeat(74));

    let mut rng = Rng(0xDEE9);
    let mut shallow_stronger = 0usize;
    for name in BANDS {
        let (so, sp, sx) = analyse(name, &sd, &sl, &mut rng);
        let (do_, dp, dx) = analyse(name, &dd, &dl, &mut rng);
        let period = doodson::constituent(name).unwrap().period_days();
        println!(
            "{:<6} {:>8.3} | {:>8.2}% {:>7.2}% {:>6.4}{} | {:>8.2}% {:>7.2}% {:>6.4}{}",
            name,
            period,
            eps(so, sd.len()) * 100.0,
            eps(sx, sd.len()) * 100.0,
            sp,
            if sp < 0.05 { "*" } else { " " },
            eps(do_, dd.len()) * 100.0,
            eps(dx, dd.len()) * 100.0,
            dp,
            if dp < 0.05 { "*" } else { " " },
        );
        if eps(so, sd.len()) > eps(do_, dd.len()) {
            shallow_stronger += 1;
        }
    }

    println!(
        "\nshallow exceeded deep in {shallow_stronger}/{} bands (chance: {:.1})",
        BANDS.len(),
        BANDS.len() as f64 / 2.0
    );
    println!(
        "eps = observed fractional modulation; lim = what we could have detected.\n\
         A sign test on {}/{} is p = {:.3} one-tailed -- reported for completeness,\n\
         not as a result, since the bands are not independent of each other.",
        shallow_stronger,
        BANDS.len(),
        // Binomial tail P(X >= k) for p=0.5, n=6.
        {
            let n = BANDS.len() as u32;
            let k = shallow_stronger as u32;
            let c = |n: u32, r: u32| -> f64 {
                (0..r).fold(1.0, |a, i| a * (n - i) as f64 / (i + 1) as f64)
            };
            (k..=n).map(|i| c(n, i)).sum::<f64>() / 2f64.powi(n as i32)
        }
    );
}
