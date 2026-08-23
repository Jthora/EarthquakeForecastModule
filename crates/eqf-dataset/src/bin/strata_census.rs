//! How much independent data does each stratum actually hold?
//!
//!     cargo run --release -p eqf-dataset --bin strata_census
//!
//! Run before pre-registering the stratified analysis, because the answer decides
//! whether that analysis is worth pre-registering. Splitting the sample raises the
//! detectable effect as 1/sqrt(n), and the global scan already needed 10,813
//! independent events to reach 0.04 log-odds per SD. A stratum with a tenth of
//! that reaches only 0.13, which is far above any effect the tidal-triggering
//! literature reports -- and a test that cannot detect what it is looking for adds
//! a multiplicity penalty while contributing nothing.
//!
//! Counting events is not testing them: no feature is touched here.

use eqf::{comcat, gcmt};
use eqf_dataset::{cells::Grid, decluster, strata::{self, Depth, Mechanism}};

/// Global scan reference point: 10,813 independent strata gave a detection
/// threshold of 0.04 log-odds per SD (docs/23-ml-results.md section 7).
const REF_N: f64 = 10_813.0;
const REF_BETA: f64 = 0.04;

fn thin(events: &[gcmt::Cmt], km: f64, days: f64) -> Vec<gcmt::Cmt> {
    let mut order: Vec<usize> = (0..events.len()).collect();
    order.sort_by(|&a, &b| events[a].day.partial_cmp(&events[b].day).unwrap());
    let mut kept: Vec<usize> = Vec::new();
    for &i in &order {
        let e = &events[i];
        let mut ok = true;
        for &j in kept.iter().rev() {
            let k = &events[j];
            if e.day - k.day > days {
                break;
            }
            if decluster::distance_km(e.lat_deg, e.lon_deg, k.lat_deg, k.lon_deg) <= km {
                ok = false;
                break;
            }
        }
        if ok {
            kept.push(i);
        }
    }
    kept.into_iter().map(|i| events[i]).collect()
}

fn main() {
    let ndk = std::fs::read_to_string("data/gcmt/gcmt.ndk").expect("gcmt.ndk");
    let all = gcmt::parse_ndk(&ndk);
    println!("GCMT: {} solutions", all.len());

    // GCMT is complete globally at about M5.5 from 1976. Below that the catalogue
    // grows with the network -- but the time-stratified design compares an event
    // only against other dates in its own calendar month, which conditions out any
    // trend however large. That is what made ComCat M4.0+ usable despite a 294-fold
    // over-dispersion in annual counts. So the threshold is a power choice here,
    // not a validity one, and it is worth sweeping.
    let min_mw: f64 = std::env::args().nth(1).unwrap_or("5.5".into()).parse().unwrap();
    let sel: Vec<gcmt::Cmt> = all
        .iter()
        .copied()
        .filter(|c| c.mw >= min_mw && c.day >= -8766.0 && c.day <= 9132.0)
        .collect();
    println!("M{min_mw}+ 1976-2024: {} solutions", sel.len());

    // Decluster on the ComCat-equivalent representation.
    let as_quakes: Vec<comcat::Quake> = sel
        .iter()
        .map(|c| comcat::Quake {
            day: c.day,
            lat_deg: c.lat_deg,
            lon_deg: c.lon_deg,
            depth_km: c.depth_km,
            magnitude: c.mw,
        })
        .collect();
    let (main, _) = decluster::gardner_knopoff(&as_quakes);
    let keep: std::collections::HashSet<u64> =
        main.iter().map(|q| q.day.to_bits()).collect();
    let indep: Vec<gcmt::Cmt> = sel
        .iter()
        .copied()
        .filter(|c| keep.contains(&c.day.to_bits()))
        .collect();
    println!("declustered: {}", indep.len());

    let grid = Grid::new(100.0);
    let _ = &grid;

    println!(
        "\n{:<26} {:>9} {:>9} {:>9}   {}",
        "stratum", "declust", "thinned", "beta_min", "verdict"
    );
    let report = |name: &str, set: &[gcmt::Cmt]| {
        let t = thin(set, 500.0, 365.0);
        let n = t.len() as f64;
        let beta = if n > 0.0 { REF_BETA * (REF_N / n).sqrt() } else { f64::INFINITY };
        // Published tidal-triggering effects on ordinary crust sit at a few
        // percent; this programme bounded M2 below 3.88% by five routes.
        let verdict = if beta <= 0.05 {
            "can test few-percent effects"
        } else if beta <= 0.12 {
            "underpowered for a few percent"
        } else {
            "too small to test anything"
        };
        println!(
            "{name:<26} {:>9} {:>9} {:>9.3}   {verdict}",
            set.len(),
            t.len(),
            beta
        );
    };

    report("ALL (reference)", &indep);
    println!();
    for m in [Mechanism::Thrust, Mechanism::Normal, Mechanism::StrikeSlip] {
        let s: Vec<gcmt::Cmt> = indep.iter().copied().filter(|c| strata::mechanism(c) == m).collect();
        report(m.name(), &s);
    }
    println!();
    for d in [Depth::Shallow, Depth::Deep] {
        let s: Vec<gcmt::Cmt> = indep.iter().copied().filter(|c| strata::depth_class(c) == d).collect();
        report(if d == Depth::Shallow { "shallow (<70 km)" } else { "deep (>=70 km)" }, &s);
    }
    println!();
    for m in [Mechanism::Thrust, Mechanism::Normal, Mechanism::StrikeSlip] {
        let s: Vec<gcmt::Cmt> = indep
            .iter()
            .copied()
            .filter(|c| strata::mechanism(c) == m && strata::depth_class(c) == Depth::Shallow)
            .collect();
        report(&format!("shallow {}", m.name()), &s);
    }
}
