//! B1: reproduce Métivier et al. (2009) — external validation of the stress path.
//!
//!     cargo run --release --example metivier_reproduction
//!
//! # Why this specific paper
//!
//! Every check in this project has been my code against my code. The moonquake
//! test validated *timing*; the Ide size-dependence was a by-product of trying to
//! disprove ourselves. Nothing has validated the **stress** path against someone
//! else's published number.
//!
//! Métivier, de Viron, Conrad, Renault, Diament & Patau (2009), *EPSL* 278,
//! 370–375, analysed **442,412 NEIC events** and reported:
//!
//! > Earthquakes occur slightly more often at the time of **ground uplift** by the
//! > Earth tide, when normal stresses are reduced within the lithosphere.
//!
//! At ~99% confidence, with the anomaly **larger for smaller and shallower**
//! events.
//!
//! # The test
//!
//! Ground uplift tracks the tide-generating potential. For a degree-2 solid
//! harmonic `V ∝ r²`, so `∂²V/∂r² = 2V/R²` — meaning the **Down-Down component of
//! the local tidal tensor is proportional to V**, and uplift maximum is `T_DD`
//! maximum. No Love numbers or fault geometry are needed: the sign is all that
//! matters, and it is geometry-free.
//!
//! So: **do more than half of events occur while `T_DD > 0`?**
//!
//! # Making it tractable
//!
//! 488k events x 400 null trials is ~390M ephemeris calls if done naively. But the
//! tidal tensor depends only on **time**, not location — the degree-2 field is five
//! numbers globally. So tensors are computed once on a fine time grid and each
//! event only interpolates and rotates to its own latitude and longitude. That is
//! the architecture `docs/06-engine-architecture.md` §1 describes, finally earning
//! its place.
//!
//! Observed and null use the identical path, so interpolation error affects both
//! equally and the comparison stays fair.
//!
//! Catalogue is ComCat (NEIC is its source) at M ≥ 4.0 from 1976, giving a count
//! close to theirs. Their depth and magnitude dependences are checked too, since
//! our own P3.5 found **no** depth dependence and this is an independent look at
//! the same claim.

use eqf::comcat;
use ph_core::{fault, field::TidalField};
use rustspice_core::{Et, KernelSet};

const NULL_TRIALS: usize = 400;
const BLOCK_DAYS: f64 = 30.0;

struct Rng(u64);
impl Rng {
    fn next_f64(&mut self) -> f64 {
        self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        ((self.0 >> 11) as f64) / ((1u64 << 53) as f64)
    }
}

fn main() -> rustspice_core::Result<()> {
    let mut ks = KernelSet::new();
    for k in ["naif0012.tls", "de440s.bsp", "pck00011.tpc", "gm_de440.tpc"] {
        ks.add_file(format!("kernels/{k}"))?;
    }
    let mut spice = ks.open()?;

    let ev = comcat::select(
        &comcat::parse_catalog(
            &std::fs::read_to_string("data/comcat/global_m40.csv").expect("run fetch-comcat 4.0"),
        ),
        4.0,
        None,
    );
    println!(
        "ComCat M4.0+: {} events, {:.1} yr   (Metivier: 442,412 NEIC events)",
        ev.len(),
        (ev[ev.len() - 1].day - ev[0].day) / 365.25
    );

    let epoch2000 = spice.parse_time("2000-01-01T00:00:00")?;
    let field = TidalField::on_earth(&mut spice, "IAU_EARTH")?;

    // Tensor grid, computed once. 15 minutes is ~7 degrees of M2 phase, so linear
    // interpolation is well inside the precision a sign test needs.
    let step = 15.0 / 1440.0;
    let (g0, g1) = (ev[0].day - 2.0, ev[ev.len() - 1].day + 2.0);
    let ng = ((g1 - g0) / step) as usize + 1;
    println!("precomputing {ng} tensors on a {:.0}-minute grid...", step * 1440.0);
    let grid_epochs: Vec<Et> = (0..ng)
        .map(|i| Et(epoch2000.0 + (g0 + i as f64 * step) * 86400.0))
        .collect();
    let grid = field.tensors(&mut spice, &grid_epochs)?;
    println!("  done; {:.0} MB", (grid.len() * 72) as f64 / 1e6);

    // T_DD at an event's own location, by interpolating the grid then rotating.
    // Positive means uplift, since T_DD is proportional to the potential.
    let dd_at = |day: f64, lat: f64, lon: f64| -> f64 {
        let x = (day - g0) / step;
        let i = (x.floor() as usize).min(grid.len() - 2);
        let f = x - i as f64;
        let mut m = ph_core::tidal::TidalTensor::default();
        for a in 0..3 {
            for b in 0..3 {
                m.m[a][b] = grid[i].m[a][b] * (1.0 - f) + grid[i + 1].m[a][b] * f;
            }
        }
        fault::to_local_ned(&m, lat, lon).m[2][2]
    };

    let dd = |sub: &[comcat::Quake], shift: &[f64]| -> Vec<f64> {
        sub.iter()
            .zip(shift)
            .map(|(q, s)| dd_at(q.day + s, q.lat_deg, q.lon_deg))
            .collect()
    };

    let zero = vec![0.0; ev.len()];
    let obs = dd(&ev, &zero);
    let frac = obs.iter().filter(|&&x| x > 0.0).count() as f64 / obs.len() as f64;

    // Block-shift null: aftershock clustering makes the binomial 50% invalid.
    let t0 = ev[0].day;
    let nb = (((ev[ev.len() - 1].day - t0) / BLOCK_DAYS).floor() as usize) + 2;
    let mut rng = Rng(0x4E1C);
    let mut null = Vec::with_capacity(NULL_TRIALS);
    println!("running {NULL_TRIALS} null trials over {} events...", ev.len());
    for _ in 0..NULL_TRIALS {
        let offs: Vec<f64> = (0..nb).map(|_| rng.next_f64() * 1.03505).collect();
        let shift: Vec<f64> = ev
            .iter()
            .map(|q| offs[(((q.day - t0) / BLOCK_DAYS).floor() as usize).min(offs.len() - 1)])
            .collect();
        let v = dd(&ev, &shift);
        null.push(v.iter().filter(|&&x| x > 0.0).count() as f64 / v.len() as f64);
    }
    null.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let ge = null.iter().filter(|&&x| x >= frac).count();
    let p = (ge as f64 + 1.0) / (NULL_TRIALS as f64 + 1.0);

    println!("\n=== uplift preference ===");
    println!("  events with T_DD > 0 (uplift): {:.4}%", frac * 100.0);
    println!("  null median                  : {:.4}%", null[null.len() / 2] * 100.0);
    println!("  null max                     : {:.4}%", null[null.len() - 1] * 100.0);
    println!("  p = {p:.4}{}", if p < 0.05 { "  *" } else { "" });
    println!(
        "  Metivier: events occur MORE often at uplift, ~99% confidence -> {}",
        if p < 0.05 && frac > null[null.len() / 2] {
            "REPRODUCED"
        } else if frac > null[null.len() / 2] {
            "same sign, not significant here"
        } else {
            "NOT reproduced -- opposite sign"
        }
    );

    // Their secondary claims: stronger for shallower and for smaller events.
    println!("\n=== dependences they report (shallower and smaller = stronger) ===");
    for (label, key) in [("depth (km)", 0u8), ("magnitude", 1)] {
        let mut rows: Vec<(f64, f64)> = ev
            .iter()
            .zip(&obs)
            .map(|(q, &d)| (if key == 0 { q.depth_km } else { q.magnitude }, d))
            .filter(|(k, _)| k.is_finite())
            .collect();
        rows.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        let per = rows.len() / 4;
        print!("  {label:<12}");
        for q in 0..4 {
            let hi = if q == 3 { rows.len() } else { (q + 1) * per };
            let s = &rows[q * per..hi];
            let f = s.iter().filter(|(_, d)| *d > 0.0).count() as f64 / s.len() as f64;
            print!("  Q{}={:.3}%", q + 1, f * 100.0);
        }
        println!();
    }
    println!("  (our P3.5 found no depth dependence; this is an independent look)");
    Ok(())
}
