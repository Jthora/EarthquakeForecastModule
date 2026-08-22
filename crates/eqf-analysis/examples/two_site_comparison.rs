//! C5: Parkfield versus Cascadia — the independence test.
//!
//!     cargo run --release --example two_site_comparison
//!
//! Every Parkfield result rests on one location with co-located families, so its
//! internal consistency was a coherence check rather than confirmation. Cascadia
//! differs in tectonic setting (subduction thrust versus strike-slip transform),
//! in geography (~1000 km away), and — importantly — in **detection method**
//! (envelope cross-correlation versus template matching).
//!
//! That last difference is what makes the comparison worth something. Parkfield's
//! diurnal artifact (D1: S1 power 16,245 against a null expectation of 1) comes
//! from template matching against a time-varying noise floor. A different pipeline
//! should carry a *different* artifact, so agreement between the sites is hard to
//! explain instrumentally.
//!
//! Phase-only, so no fault geometry is needed and the two are directly comparable.

use eqf::{cascadia, parkfield};
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

fn power(c: &doodson::Constituent, times: &[f64], shift: impl Fn(f64) -> f64) -> f64 {
    let (mut a, mut b) = (0.0f64, 0.0f64);
    for &t in times {
        let (s, co) = c.phase_at(t + shift(t)).sin_cos();
        a += co;
        b += s;
    }
    (a * a + b * b) / times.len() as f64
}

/// Observed power, null median, and empirical p for one constituent.
fn analyse(name: &str, times: &[f64], rng: &mut Rng) -> (f64, f64, f64) {
    let c = doodson::constituent(name).unwrap();
    let period = c.period_days();
    let block = (4.0 * period).max(30.0);
    let t0 = times[0];
    let t1 = times[times.len() - 1];
    let n_blocks = (((t1 - t0) / block).floor() as usize) + 2;

    let observed = power(c, times, |_| 0.0);
    let mut null: Vec<f64> = (0..NULL_TRIALS)
        .map(|_| {
            let offs: Vec<f64> = (0..n_blocks).map(|_| rng.next_f64() * period).collect();
            power(c, times, |t| {
                offs[(((t - t0) / block).floor() as usize).min(offs.len() - 1)]
            })
        })
        .collect();
    null.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let ge = null.iter().filter(|&&x| x >= observed).count();
    (
        observed,
        null[null.len() / 2],
        (ge as f64 + 1.0) / (NULL_TRIALS as f64 + 1.0),
    )
}

fn load_sorted<T, F: Fn(&T) -> f64>(v: &[T], f: F) -> Vec<f64> {
    let mut t: Vec<f64> = v.iter().map(&f).collect();
    t.sort_by(|a, b| a.partial_cmp(b).unwrap());
    t
}

fn main() {
    let pk = parkfield::parse_catalog(
        &std::fs::read_to_string("data/parkfield/LFEcat_Apr2001-Apr2024.csv")
            .expect("run scripts/fetch-parkfield.sh"),
    );
    let cs = cascadia::parse_catalog(
        &std::fs::read_to_string("data/cascadia/cascadia_tremor.csv")
            .expect("run scripts/fetch-cascadia.sh"),
    );

    let pk_t = load_sorted(&pk, |e| e.day);
    let cs_t = load_sorted(&cs, |e| e.day);
    println!(
        "Parkfield LFEs  : {:>9} events, {:.1} yr",
        pk_t.len(),
        (pk_t[pk_t.len() - 1] - pk_t[0]) / 365.25
    );
    println!(
        "Cascadia tremor : {:>9} events, {:.1} yr\n",
        cs_t.len(),
        (cs_t[cs_t.len() - 1] - cs_t[0]) / 365.25
    );

    println!(
        "{:<6} {:>9} | {:>10} {:>9} {:>8} | {:>10} {:>9} {:>8}",
        "band", "period", "PK D2/N", "PK ratio", "PK p", "CS D2/N", "CS ratio", "CS p"
    );
    println!("{}", "-".repeat(82));

    let mut rng = Rng(0xC5C0DE);
    let mut agree = 0usize;
    let mut both = 0usize;
    for name in BANDS {
        let (po, pm, pp) = analyse(name, &pk_t, &mut rng);
        let (co, cm, cp) = analyse(name, &cs_t, &mut rng);
        let period = doodson::constituent(name).unwrap().period_days();
        println!(
            "{:<6} {:>9.4} | {:>10.1} {:>9.1} {:>8.4}{} | {:>10.1} {:>9.1} {:>8.4}{}",
            name,
            period,
            po,
            po / pm,
            pp,
            if pp < 0.05 { "*" } else { " " },
            co,
            co / cm,
            cp,
            if cp < 0.05 { "*" } else { " " }
        );
        if (pp < 0.05) == (cp < 0.05) {
            agree += 1;
        }
        if pp < 0.05 && cp < 0.05 {
            both += 1;
        }
    }

    println!(
        "\nagreement: {agree}/{} constituents give the same verdict; {both} significant at both sites",
        BANDS.len()
    );
    println!("null floor 1/(n+1) = {:.4}", 1.0 / (NULL_TRIALS as f64 + 1.0));
    println!(
        "\nDetection differs between the sites (template matching vs envelope\n\
         cross-correlation), so a shared artifact is an unlikely explanation for\n\
         any constituent significant at both."
    );
}
