//! Catalogue inventory: what survives declustering, and what the matched design costs.
//!
//!     cargo run --release -p eqf-dataset --bin inspect -- data/comcat/global_m40.csv
//!
//! Prints the numbers that decide the experiment's shape — how many independent
//! events there are to learn from, and how large the design matrix would be.

use eqf::comcat;
use eqf_dataset::{cells::Grid, decluster, sampling};

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let path = args.first().map(String::as_str).unwrap_or("data/comcat/global_m40.csv");
    let csv = std::fs::read_to_string(path).expect("read catalogue");
    let all = comcat::parse_catalog(&csv);
    println!("catalogue {path}\n  parsed {} events", all.len());

    let (lo, hi) = all.iter().fold((f64::MAX, f64::MIN), |(a, b), q| (a.min(q.day), b.max(q.day)));
    println!("  span   day {lo:.1} to {hi:.1}  ({:.1} years)", (hi - lo) / 365.25);

    let grid = Grid::new(100.0);
    println!("  grid   {} cells of ~100 km", grid.len());

    for min_mag in [4.0, 4.5, 5.0, 5.5, 6.0] {
        let sel: Vec<_> = all.iter().copied().filter(|q| q.magnitude >= min_mag).collect();
        if sel.is_empty() {
            continue;
        }
        let t0 = std::time::Instant::now();
        let (main, dep) = decluster::gardner_knopoff(&sel);
        let secs = t0.elapsed().as_secs_f64();

        let occupied: std::collections::HashSet<usize> =
            main.iter().map(|q| grid.cell(q.lat_deg, q.lon_deg)).collect();

        // Design-matrix size for the full feature vector at 10 controls per case.
        const FEATURES: usize = 9_816;
        let rows = main.len() * 11;
        let gb = rows as f64 * FEATURES as f64 * 4.0 / 1e9;

        println!(
            "\n  M{min_mag}+   {:>7} events -> {:>7} independent ({:.1}% removed) in {secs:.1}s",
            sel.len(),
            main.len(),
            100.0 * dep.len() as f64 / sel.len() as f64
        );
        println!(
            "           {:>7} cells occupied, {:.1} events per occupied cell",
            occupied.len(),
            main.len() as f64 / occupied.len() as f64
        );
        println!("           design matrix {rows} rows x {FEATURES} f32 = {gb:.1} GB");
    }

    // Build the primary design once, on the declustered M5.5 set, and report the
    // matched-set statistics that determine how much the model has to work with.
    let sel: Vec<_> = all.iter().copied().filter(|q| q.magnitude >= 5.5).collect();
    let (main, _) = decluster::gardner_knopoff(&sel);
    let mut rng = sampling::Rng::seed(20260822);
    let rows = sampling::build(
        &main,
        &grid,
        sampling::Scheme::DayOffset { max_days: 5 },
        10,
        (lo, hi),
        &mut rng,
    );
    let cases = rows.iter().filter(|r| r.case).count();
    println!(
        "\n  primary design (M5.5+, declustered, +/-1..5 day offsets, 10 controls)"
    );
    println!("           {cases} cases, {} controls, {} rows total", rows.len() - cases, rows.len());
    let short = rows
        .chunk_by(|a, b| a.stratum == b.stratum)
        .filter(|s| s.len() < 11)
        .count();
    println!("           {short} strata got fewer than 10 controls (span edges)");
}
