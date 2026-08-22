//! A5/A6: per-nest fault-orientation search against Weber's criterion.
//!
//!     cargo run --release --example moonquake_nest_coulomb
//!
//! Pooled-catalogue phase tests cannot produce a falsifiable claim here: the
//! catalogue and the forcing share the anomalistic month, so clustering is
//! guaranteed at *some* phase (see `examples/moonquake_tidal_phase.rs`). Weber,
//! Bills & Johnson (2009) work per nest instead, and so does this.
//!
//! **Criterion.** For each nest, search fault orientations for the one where
//! ΔCFS at event times "best approximates a constant" — i.e. minimises
//!
//! ```text
//! score = std(ΔCFS at events) / std(ΔCFS over the whole span)
//! ```
//!
//! Under random timing the expectation is ~1. Lower means events recur at a
//! consistent stress state.
//!
//! **The null matters more than the score.** Minimising over thousands of planes
//! will find a low score by chance, so the comparison is against the *same search*
//! run on random event times — best-of-grid versus best-of-grid.

use eqf::{apollo};
use ph_core::{fault, field::TidalField, tidal::TidalTensor};
use rustspice_core::{Et, KernelSet};

const MIN_EVENTS: usize = 20;
// BH-FDR over 74 nests needs p <= 0.05/74 = 6.8e-4 for the strongest nest to
// survive. The empirical floor is 1/(n+1), so anything under ~1470 trials cannot
// clear it no matter how strong the signal.
const NULL_TRIALS: usize = 2000;
const MU_GRID: &[f64] = &[0.0, 0.2, 0.4, 0.6, 0.8];

/// Deterministic LCG — reproducible nulls without a dependency.
struct Rng(u64);
impl Rng {
    fn next_usize(&mut self, n: usize) -> usize {
        self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        ((self.0 >> 33) as usize) % n
    }
}

fn main() -> rustspice_core::Result<()> {
    let mut ks = KernelSet::new();
    for k in [
        "naif0012.tls", "de440s.bsp", "pck00011.tpc",
        "gm_de440.tpc", "moon_pa_de440_200625.bpc", "moon_de440_250416.tf",
    ] {
        ks.add_file(format!("kernels/{k}"))?;
    }
    let mut spice = ks.open()?;

    let events = apollo::parse_levent(
        &std::fs::read_to_string("data/apollo/levent.1008weber.csv").expect("catalogue"),
    );
    let locations = apollo::parse_dm_locations(
        &std::fs::read_to_string("data/apollo/nakamura_2005_dm_locations.csv").expect("locations"),
    );
    println!("nest locations: {}", locations.len());

    let all = apollo::deep_moonquakes(&events, None);
    let days = all.times();
    let (t0, t1) = (days[0], days[days.len() - 1]);

    // Reference grid, computed once in MOON_PA and reused for every nest.
    let step = 0.25;
    let n = ((t1 - t0) / step) as usize;
    let grid_days: Vec<f64> = (0..n).map(|i| t0 + i as f64 * step).collect();
    let epoch0 = spice.parse_time("1969-01-01T00:00:00")?;
    let epochs: Vec<Et> = grid_days.iter().map(|&d| Et(epoch0.0 + d * 86400.0)).collect();

    let lunar = TidalField::on_moon(&mut spice)?;
    let reference = lunar.tensors(&mut spice, &epochs)?;
    println!("reference grid: {} samples at {step} d\n", reference.len());

    let planes = fault::plane_grid(15.0, 15.0, 30.0);
    println!("search space: {} planes x {} friction values\n", planes.len(), MU_GRID.len());

    // Coefficients depend only on geometry and friction, so build them once.
    // They involve trig; recomputing inside the null loop would dominate runtime.
    let coeffs: Vec<[f64; 6]> = planes
        .iter()
        .flat_map(|p| MU_GRID.iter().map(move |&mu| fault::coulomb_coefficients(p, mu)))
        .collect();

    // Index of the nearest reference sample to a given day.
    let idx_of = |d: f64| (((d - t0) / step).round() as usize).min(reference.len() - 1);

    println!(
        "{:>5} {:>6} {:>9} {:>9} {:>8} {:>9} {:>8}",
        "nest", "events", "score", "unif_best", "p_unif", "shift_best", "p_shift"
    );
    println!("{}", "-".repeat(64));

    let mut rows = Vec::new();
    for loc in &locations {
        let nest_days = apollo::deep_moonquakes(&events, Some(loc.nest)).times();
        if nest_days.len() < MIN_EVENTS {
            continue;
        }

        // Rotate the shared reference into this nest's local frame.
        let local: Vec<TidalTensor> = reference
            .iter()
            .map(|t| fault::to_local_ned(t, loc.lat_deg, loc.lon_deg))
            .collect();
        let cov_ref = fault::component_covariance(&local);
        // Denominators are fixed for this nest; precompute alongside the coeffs.
        let denoms: Vec<f64> = coeffs.iter().map(|c| fault::coulomb_std(&cov_ref, c)).collect();

        let ev_idx: Vec<usize> = nest_days.iter().map(|&d| idx_of(d)).collect();
        let ev_tensors: Vec<TidalTensor> = ev_idx.iter().map(|&i| local[i]).collect();
        let cov_ev = fault::component_covariance(&ev_tensors);

        // Best score over the whole search space.
        let best = |cov: &[[f64; 6]; 6]| -> f64 {
            let mut b = f64::MAX;
            for (c, &denom) in coeffs.iter().zip(&denoms) {
                if denom <= 0.0 {
                    continue;
                }
                let s = fault::coulomb_std(cov, c) / denom;
                if s < b {
                    b = s;
                }
            }
            b
        };

        let observed = best(&cov_ev);

        // Two nulls, deliberately different in what they preserve.
        //
        // UNIFORM draws random times, destroying the catalogue's temporal
        // clustering. That makes it test "are events clustered in time at all?",
        // which is trivially yes -- events at nearby times see similar tensors, so
        // their variance is low regardless of any tidal alignment.
        //
        // SHIFT slides the whole nest sequence by a random offset, preserving
        // relative spacing and breaking only the alignment with the forcing. That
        // is the question we actually mean to ask.
        let mut rng = Rng(0x5EED ^ (loc.nest as u64) << 8);
        let mut null = Vec::with_capacity(NULL_TRIALS);
        let mut null_shift = Vec::with_capacity(NULL_TRIALS);
        for _ in 0..NULL_TRIALS {
            let sample: Vec<TidalTensor> = (0..nest_days.len())
                .map(|_| local[rng.next_usize(local.len())])
                .collect();
            null.push(best(&fault::component_covariance(&sample)));

            let off = rng.next_usize(local.len());
            let shifted: Vec<TidalTensor> = ev_idx
                .iter()
                .map(|&i| local[(i + off) % local.len()])
                .collect();
            null_shift.push(best(&fault::component_covariance(&shifted)));
        }
        null_shift.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let better_shift = null_shift.iter().filter(|&&x| x <= observed).count();
        let p_shift = (better_shift as f64 + 1.0) / (NULL_TRIALS as f64 + 1.0);
        null.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let better = null.iter().filter(|&&x| x <= observed).count();
        let p = (better as f64 + 1.0) / (NULL_TRIALS as f64 + 1.0);

        println!(
            "{:>5} {:>6} {:>9.4} {:>9.4} {:>8.4} {:>9.4} {:>8.4}{}",
            loc.nest,
            nest_days.len(),
            observed,
            null[0],
            p,
            null_shift[0],
            p_shift,
            if p_shift < 0.05 { "  *" } else { "" }
        );
        rows.push((loc.nest, p, p_shift));
    }

    let sig_u = rows.iter().filter(|(_, p, _)| *p < 0.05).count();
    let sig_s = rows.iter().filter(|(_, _, p)| *p < 0.05).count();
    println!(
        "\nuniform null: {sig_u}/{} nests below p=0.05",
        rows.len()
    );
    println!(
        "shift null:   {sig_s}/{} nests below p=0.05   (expected by chance: {:.1})",
        rows.len(),
        0.05 * rows.len() as f64
    );
    println!("null floor is 1/(n+1) = {:.4}", 1.0 / (NULL_TRIALS as f64 + 1.0));

    // Benjamini-Hochberg. 74 nests tested, so an uncorrected count is not a
    // result on its own.
    let mut ps: Vec<f64> = rows.iter().map(|(_, _, p)| *p).collect();
    ps.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let m = ps.len() as f64;
    let mut k_max = 0usize;
    for (i, &p) in ps.iter().enumerate() {
        if p <= 0.05 * (i as f64 + 1.0) / m {
            k_max = i + 1;
        }
    }
    println!("\nBenjamini-Hochberg at FDR 0.05: {k_max}/{} nests survive", ps.len());

    // Ensemble test: is the count of nominally significant nests itself unusual?
    let obs = rows.iter().filter(|(_, _, p)| *p < 0.05).count() as f64;
    let mean = 0.05 * m;
    let sd = (m * 0.05 * 0.95).sqrt();
    println!(
        "ensemble: {obs:.0} nominal hits vs {mean:.1} expected, sd {sd:.2} -> {:.1} sigma",
        (obs - mean) / sd
    );
    Ok(())
}
