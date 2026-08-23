//! Declustered event counts per calendar year, for the Poisson dispersion check.

use eqf::comcat;
use eqf_dataset::{decluster, sampling::civil_from_days};

fn main() {
    let min_mag: f64 = std::env::args().nth(1).unwrap_or("5.5".into()).parse().unwrap();
    let csv = std::fs::read_to_string("data/comcat/global_m40.csv").expect("catalogue");
    let all = comcat::parse_catalog(&csv);
    let sel: Vec<_> = all.into_iter().filter(|q| q.magnitude >= min_mag).collect();
    let (main, _) = decluster::gardner_knopoff(&sel);
    let mut counts: std::collections::BTreeMap<i64, usize> = Default::default();
    for q in &main {
        let (y, _, _) = civil_from_days(q.day.floor() as i64);
        *counts.entry(y).or_insert(0) += 1;
    }
    for (y, c) in counts {
        println!("{y} {c}");
    }
}
