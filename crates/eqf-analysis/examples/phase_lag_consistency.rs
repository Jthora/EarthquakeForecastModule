//! H1/H2: is the preferred phase consistent with a single physical mechanism?
//!
//!     cargo run --release --example phase_lag_consistency
//!
//! # What this tests
//!
//! If triggering is real, the preferred phase at each constituent should reflect
//! **one** physical relationship, not several unrelated numbers. Measuring the
//! observed preferred phase *relative to the ΔCFS phase at that site* gives a lag
//! per constituent, and its behaviour across constituents is diagnostic:
//!
//! | Lag pattern | Mechanism |
//! |---|---|
//! | ≈ 0° at every constituent | **amplitude-driven** — failure at peak ΔCFS |
//! | ≈ ±90° at every constituent | **rate-driven** — failure at peak dΔCFS/dt |
//! | scales as 1/period | a **constant time delay** — nucleation |
//! | scattered | no single mechanism; something is wrong |
//!
//! **This subsumes H2.** For a single constituent `dΔCFS/dt` is just ΔCFS shifted
//! 90°, and Schuster power is rotation-invariant, so testing "rate versus
//! amplitude" by power is vacuous — the two are identical. The *phase* is where
//! the distinction lives, so a lag near ±90° is the rate-driven signature.
//!
//! Only constituents that survived C4 at both sites are used (M2, N2, O1),
//! plus Q1 which survived at Parkfield alone. Long-period constituents are
//! excluded: their phases are not significant, so their lags are noise.

use eqf::parkfield;
use ph_core::{doodson, fault, field::TidalField, love::Elastic};
use rustspice_core::{Et, KernelSet};

const SAF: fault::FaultPlane = fault::FaultPlane {
    strike_deg: 137.0,
    dip_deg: 90.0,
    rake_deg: 180.0,
};
const BANDS: &[&str] = &["M2", "N2", "O1", "Q1"];
const STEP: f64 = 0.02;

fn norm180(deg: f64) -> f64 {
    let mut d = deg % 360.0;
    if d > 180.0 {
        d -= 360.0;
    }
    if d < -180.0 {
        d += 360.0;
    }
    d
}

fn main() -> rustspice_core::Result<()> {
    let mut ks = KernelSet::new();
    for k in ["naif0012.tls", "de440s.bsp", "pck00011.tpc", "gm_de440.tpc"] {
        ks.add_file(format!("kernels/{k}"))?;
    }
    let mut spice = ks.open()?;

    let ev = parkfield::parse_catalog(
        &std::fs::read_to_string("data/parkfield/LFEcat_Apr2001-Apr2024.csv").expect("parkfield"),
    );
    let mut times: Vec<f64> = ev.iter().map(|e| e.day).collect();
    times.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let (t0, t1) = (times[0], times[times.len() - 1]);
    let (lat, lon, _) =
        parkfield::family_location(&ev, &parkfield::families(&ev)[0].0).unwrap();
    println!("Parkfield: {} events, {:.1} yr\n", times.len(), (t1 - t0) / 365.25);

    // Solid-tide dCFS series at the site.
    let n = ((t1 - t0) / STEP) as usize;
    let epoch2000 = spice.parse_time("2000-01-01T00:00:00")?;
    let days: Vec<f64> = (0..n).map(|i| t0 + i as f64 * STEP).collect();
    let epochs: Vec<Et> = days.iter().map(|&d| Et(epoch2000.0 + d * 86400.0)).collect();
    let tensors = TidalField::on_earth(&mut spice, "IAU_EARTH")?.tensors(&mut spice, &epochs)?;
    let e = Elastic::EARTH;
    let cfs: Vec<f64> = tensors
        .iter()
        .map(|t| e.stress(fault::coulomb(&fault::to_local_ned(t, lat, lon), &SAF, 0.4)))
        .collect();

    println!(
        "{:<5} {:>9} {:>11} {:>11} {:>10} {:>12}",
        "band", "period_h", "phase_cfs", "phase_obs", "lag_deg", "lag_minutes"
    );
    println!("{}", "-".repeat(62));

    let mut rows = Vec::new();
    for name in BANDS {
        let c = doodson::constituent(name).unwrap();
        let period = c.period_days();

        // dCFS phase at this site, by least squares on the analytic argument.
        let (mut sa, mut sb) = (0.0f64, 0.0f64);
        for (k, &d) in days.iter().enumerate() {
            let (s, co) = c.phase_at_longitude(d, lon).sin_cos();
            sa += cfs[k] * co;
            sb += cfs[k] * s;
        }
        let phase_cfs = (-sb).atan2(sa).to_degrees();

        // Observed preferred phase in the same argument.
        let (mut a, mut b) = (0.0f64, 0.0f64);
        for &t in &times {
            let (s, co) = c.phase_at_longitude(t, lon).sin_cos();
            a += co;
            b += s;
        }
        let phase_obs = b.atan2(a).to_degrees();

        let lag = norm180(phase_obs - phase_cfs);
        let minutes = lag / 360.0 * period * 24.0 * 60.0;
        println!(
            "{:<5} {:>9.2} {:>11.1} {:>11.1} {:>10.1} {:>12.1}",
            name,
            period * 24.0,
            norm180(phase_cfs),
            norm180(phase_obs),
            lag,
            minutes
        );
        rows.push((name.to_string(), period, lag, minutes));
    }

    // Group by species: semidiurnal constituents have Doodson n1 = 2, diurnal
    // n1 = 1. A pooled spread hides the structure -- the question is whether lags
    // are consistent WITHIN a band and how the bands differ.
    let spread = |v: &[f64]| {
        if v.len() < 2 {
            return 0.0;
        }
        let m = v.iter().sum::<f64>() / v.len() as f64;
        (v.iter().map(|x| (x - m).powi(2)).sum::<f64>() / v.len() as f64).sqrt()
    };
    let mean = |v: &[f64]| v.iter().sum::<f64>() / v.len() as f64;

    let semi: Vec<f64> = rows.iter().filter(|r| r.1 < 0.75).map(|r| r.2).collect();
    let diur: Vec<f64> = rows.iter().filter(|r| r.1 >= 0.75).map(|r| r.2).collect();

    println!("\nsemidiurnal lag: mean {:+.1} deg, sd {:.1}", mean(&semi), spread(&semi));
    println!("diurnal     lag: mean {:+.1} deg, sd {:.1}", mean(&diur), spread(&diur));
    println!("band difference: {:+.1} deg", mean(&diur) - mean(&semi));

    let within = spread(&semi).max(spread(&diur));
    let between = (mean(&diur) - mean(&semi)).abs();
    println!(
        "\nwithin-band spread {:.1} deg vs between-band difference {:.1} deg",
        within, between
    );
    println!(
        "\nreading: {}",
        if within > 45.0 {
            "lags scatter even within a band -- no coherent mechanism"
        } else if between < 30.0 {
            "one lag across all bands -- a single frequency-independent mechanism"
        } else {
            "lags are CONSISTENT WITHIN each band but DIFFER between them -- \n                      a frequency-dependent phase response, or a band-dependent error in the \n                      dCFS reference. Ocean loading is exactly such an error and is omitted here."
        }
    );
    Ok(())
}
