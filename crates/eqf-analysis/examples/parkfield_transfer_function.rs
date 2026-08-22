//! C4: frequency-resolved response at Parkfield — the first transfer function.
//!
//!     cargo run --release --example parkfield_transfer_function
//!
//! Every result so far sits at M2 (12.42 h). A transfer function is a *curve*, so
//! this measures response in several constituent bands separately, using complex
//! demodulation ([`ph_core::demod`]) rather than period folding.
//!
//! # Reading out a response fraction
//!
//! If rate is modulated as `1 + ε cos θ`, each event contributes mean `ε/2` to the
//! resultant, so `D²/N ≈ N ε²/4` and
//!
//! ```text
//! ε = 2 √( (D²/N) / N )
//! ```
//!
//! Dividing by the band's forcing amplitude gives the transfer function
//! `R(ω) = ε(ω) / ΔCFS_amplitude(ω)` — response per unit stress, which is the
//! quantity `docs/08` §12 is built around.
//!
//! # Exclusions, decided in advance
//!
//! **K1 (23.93 h) and S1 (24.00 h) are excluded.** D1 measured the diurnal
//! detection artifact at S1 power 16,245 against a null expectation of 1, and K1
//! sits at 1.16× that. They are degenerate with the artifact and cannot be read.
//!
//! **S2 is excluded** for the same reason: it is exactly 12.000 h and carries the
//! solar thermal tide.
//!
//! # Null: sham frequencies, not time shifts
//!
//! A first version of this used the whole-day-shift null from C2. **That was
//! wrong, and it is trap #1 recurring.** Demodulation isolates a single
//! constituent, which makes the band a near-pure tone — and against a pure tone a
//! time shift merely *rotates* the phase cluster without diluting it, so `D²` is
//! near-invariant and the null has no power. `ph_core::stats` documents exactly
//! this; the demodulation step recreated the degenerate case.
//!
//! The null used here instead runs the **identical procedure at frequencies where
//! no tidal constituent exists**. Those quiet bands share the catalogue's
//! structure, the demodulation window behaviour, and the sample size — they differ
//! only in carrying no tide. If a real constituent does not exceed them, there is
//! nothing to report.

use eqf::{parkfield};
use ph_core::{demod, fault, field::TidalField, stats};
use rustspice_core::{Et, KernelSet};

const SAF: fault::FaultPlane = fault::FaultPlane {
    strike_deg: 137.0,
    dip_deg: 90.0,
    rake_deg: 180.0,
};
const MU: f64 = 0.4;
const STEP_DAYS: f64 = 0.02;
/// Periods with no significant tidal constituent — the noise floor for the same
/// procedure. Chosen to bracket the real bands in period without colliding with
/// any named constituent.
const SHAM_PERIODS: &[(f64, f64)] = &[
    (0.72, 90.0),
    (0.83, 90.0),
    (2.6, 200.0),
    (4.1, 200.0),
    (7.3, 300.0),
    (19.4, 500.0),
    (41.0, 500.0),
    (63.0, 800.0),
    (97.0, 1000.0),
    (240.0, 1400.0),
];

/// (name, period days, demodulation window days, nearest neighbour for the beat)
const BANDS: &[(&str, f64, f64, f64)] = &[
    ("M2  semidiurnal", 0.5175, 90.0, 0.5000),
    ("O1  diurnal", 1.0758, 90.0, 0.9973),
    ("Mf  fortnightly", 13.661, 500.0, 14.765),
    ("Msf synodic fortnightly", 14.765, 500.0, 13.661),
    ("Mm  monthly", 27.555, 400.0, 14.765),
    ("Ssa semiannual", 182.62, 1400.0, 365.26),
    ("Sa  annual", 365.26, 1400.0, 182.62),
];

fn main() -> rustspice_core::Result<()> {
    let mut ks = KernelSet::new();
    for k in ["naif0012.tls", "de440s.bsp", "pck00011.tpc", "gm_de440.tpc"] {
        ks.add_file(format!("kernels/{k}"))?;
    }
    let mut spice = ks.open()?;

    let events = parkfield::parse_catalog(
        &std::fs::read_to_string("data/parkfield/LFEcat_Apr2001-Apr2024.csv")
            .expect("run scripts/fetch-parkfield.sh"),
    );
    let fams = parkfield::families(&events);
    let (lat, lon, _) = parkfield::family_location(&events, &fams[0].0).unwrap();

    let (mut lo, mut hi) = (f64::MAX, f64::MIN);
    for e in &events {
        lo = lo.min(e.day);
        hi = hi.max(e.day);
    }
    let (t0, t1) = (lo - 4500.0, hi + 4500.0);
    let n = ((t1 - t0) / STEP_DAYS) as usize;
    let epoch2000 = spice.parse_time("2000-01-01T00:00:00")?;
    let days: Vec<f64> = (0..n).map(|i| t0 + i as f64 * STEP_DAYS).collect();
    let epochs: Vec<Et> = days.iter().map(|&d| Et(epoch2000.0 + d * 86400.0)).collect();

    println!("sampling dCFS at {n} points...");
    let earth = TidalField::on_earth(&mut spice, "IAU_EARTH")?;
    let tensors = earth.tensors(&mut spice, &epochs)?;
    let cfs: Vec<f64> = tensors
        .iter()
        .map(|t| fault::coulomb(&fault::to_local_ned(t, lat, lon), &SAF, MU))
        .collect();

    let mut times: Vec<f64> = events.iter().map(|e| e.day).collect();
    times.sort_by(|a, b| a.partial_cmp(b).unwrap());
    println!("{} events\n", times.len());

    // Noise floor first: same procedure, frequencies with no tide.
    let mut sham: Vec<f64> = Vec::new();
    for (period, window) in SHAM_PERIODS {
        if let Some(b) = demod::demodulate(&days, &cfs, *period, *window) {
            let (ph, _) = b.phases(&times);
            if ph.len() >= 1000 {
                sham.push(stats::schuster(&ph, 1)[0].d_squared / ph.len() as f64);
            }
        }
    }
    sham.sort_by(|a, b| a.partial_cmp(b).unwrap());
    println!(
        "sham-frequency floor ({} quiet bands): median {:.0}, max {:.0}\n",
        sham.len(),
        sham[sham.len() / 2],
        sham.last().copied().unwrap_or(0.0)
    );

    println!(
        "{:<26} {:>8} {:>7} {:>10} {:>10} {:>11} {:>9}",
        "band", "period", "beat", "force", "D2/N", "vs sham", "p"
    );
    println!("{}", "-".repeat(86));
    for (name, period, window, neighbour) in BANDS {
        let Some(band) = demod::demodulate(&days, &cfs, *period, *window) else {
            println!("{name:<26}  demodulation failed");
            continue;
        };
        let beat = band.leakage_note(*neighbour);
        if beat > *window {
            println!(
                "{name:<26}  SKIPPED: beat {beat:.0} d exceeds window {window:.0} d"
            );
            continue;
        }

        let (phases, dropped) = band.phases(&times);
        if phases.len() < 1000 {
            println!("{name:<26}  too few events in span ({dropped} dropped)");
            continue;
        }
        let d2n = stats::schuster(&phases, 1)[0].d_squared / phases.len() as f64;
        let force = band.mean_amplitude();
        let ge = sham.iter().filter(|&&x| x >= d2n).count();
        let p = (ge as f64 + 1.0) / (sham.len() as f64 + 1.0);

        println!(
            "{:<26} {:>8.4} {:>7.0} {:>10.3e} {:>10.0} {:>11.2} {:>9.4}{}",
            name,
            period,
            beat,
            force,
            d2n,
            d2n / sham[sham.len() / 2],
            p,
            if p < 0.05 { "  *" } else { "" }
        );
    }

    println!("\nforce   = mean dCFS band amplitude (stress SHAPE, no Love numbers)");
    println!("D2/N    = Schuster power; null expectation 1 for random phases");
    println!("vs sham = D2/N relative to the median of quiet, tide-free frequencies");
    println!("p       = fraction of sham bands reaching this D2/N; floor 1/(n+1)");
    println!("\nK1, S1 and S2 excluded: degenerate with the diurnal/thermal artifact (D1)");
    println!("Time-shift nulls are NOT used here: demodulation makes each band a");
    println!("near-pure tone, against which a shift rotates phase without diluting it.");
    Ok(())
}
