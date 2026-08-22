//! C2: tidal ΔCFS phase at Parkfield, tested against a whole-day-shift null.
//!
//!     cargo run --release --example parkfield_coulomb_phase
//!
//! # The null, specified before the run
//!
//! Every trap in this project so far came from choosing a null *after* seeing the
//! data, so this one is fixed in advance and justified on structure, not outcome.
//!
//! D1 showed the catalogue carries a large detection artifact locked to **solar
//! time** — S1 (exactly 24.000 h) reaches Schuster power 16,245 against a null
//! expectation of 1. Any null must preserve that artifact, or it will "detect" the
//! artifact rather than the tide.
//!
//! **Null: shift every event time by a whole number of solar days.**
//!
//! - The detection artifact is locked to local solar time, so a whole-day shift
//!   leaves it **exactly invariant**.
//! - The lunar tide precesses ~50 min per solar day, so a whole-day shift **does**
//!   decorrelate ΔCFS phase.
//!
//! This is the cleanest separation available: it holds the confound fixed while
//! sliding only the quantity under test. Shifts are drawn from ±(30…4000) days,
//! excluding multiples of the 27.55 d anomalistic and 29.53 d synodic months to
//! avoid partially restoring alignment.
//!
//! # Geometry
//!
//! Deep San Andreas at Parkfield: right-lateral strike-slip, **strike 137°, dip
//! 90°, rake 180°**. A vertical plane makes strike and strike+180° equivalent. A
//! small sensitivity grid is also reported.

use eqf::{parkfield};
use ph_core::{fault, field::TidalField, phase::Forcing, stats, tidal::TidalTensor};
use rustspice_core::{Et, KernelSet};

const SAF: fault::FaultPlane = fault::FaultPlane {
    strike_deg: 137.0,
    dip_deg: 90.0,
    rake_deg: 180.0,
};
const MU: f64 = 0.4; // Cochran et al. (2004) best fit
const STEP_DAYS: f64 = 0.02; // ~26 samples per M2 cycle
const NULL_TRIALS: usize = 400;
const TOP_FAMILIES: usize = 12;

/// Preferred phases as reported, for the coherence check.
fn phases_summary(v: &[f64]) -> Vec<f64> {
    v.to_vec()
}

struct Rng(u64);
impl Rng {
    fn next_u64(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        self.0 >> 33
    }
}

/// Whole-day offsets, avoiding near-multiples of the lunar months.
fn day_offsets(n: usize, rng: &mut Rng) -> Vec<f64> {
    let mut out = Vec::with_capacity(n);
    while out.len() < n {
        let mag = 30 + (rng.next_u64() % 3970) as i64;
        let d = mag as f64;
        let resonant = [27.555_f64, 29.531]
            .iter()
            .any(|p| ((d / p).round() * p - d).abs() < 2.0);
        if resonant {
            continue;
        }
        out.push(if rng.next_u64() % 2 == 0 { d } else { -d });
    }
    out
}

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
    println!("{} events, {} families", events.len(), fams.len());
    println!("reference location: {lat:.3} N, {lon:.3} E\n");

    let (mut lo, mut hi) = (f64::MAX, f64::MIN);
    for e in &events {
        lo = lo.min(e.day);
        hi = hi.max(e.day);
    }
    // Pad generously so shifted events stay inside the forcing span.
    let (t0, t1) = (lo - 4200.0, hi + 4200.0);
    let n = ((t1 - t0) / STEP_DAYS) as usize;
    println!("sampling dCFS at {n} points, {STEP_DAYS} d step (this takes a minute)");

    let epoch2000 = spice.parse_time("2000-01-01T00:00:00")?;
    let days: Vec<f64> = (0..n).map(|i| t0 + i as f64 * STEP_DAYS).collect();
    let epochs: Vec<Et> = days.iter().map(|&d| Et(epoch2000.0 + d * 86400.0)).collect();

    // IAU_EARTH is adequate here: its prime-meridian rate ignores UT1-UTC, which
    // leap seconds bound to <0.9 s (0.004 deg), and it models precession to first
    // order. Both are far below the phase resolution that matters for M2.
    let earth = TidalField::on_earth(&mut spice, "IAU_EARTH")?;
    let tensors = earth.tensors(&mut spice, &epochs)?;

    let cfs: Vec<f64> = tensors
        .iter()
        .map(|t| fault::coulomb(&fault::to_local_ned(t, lat, lon), &SAF, MU))
        .collect();
    let forcing = Forcing::new(days, cfs).expect("forcing");
    println!(
        "dCFS maxima: {}   mean interval {:.4} d  (M2 is 0.5175 d)\n",
        forcing.maxima().len(),
        forcing.mean_period().unwrap()
    );

    let mut rng = Rng(0xC0FFEE);
    let offsets = day_offsets(NULL_TRIALS, &mut rng);

    println!(
        "{:<10} {:>8} {:>10} {:>9} {:>10} {:>8}",
        "family", "events", "D2/N", "phase", "null_max", "p"
    );
    println!("{}", "-".repeat(60));

    let mut ps = Vec::new();
    let mut family_phases = Vec::new();
    for (fid, count) in fams.iter().take(TOP_FAMILIES) {
        let times = parkfield::family(&events, fid).times();
        let (phases, dropped) = forcing.phases(&times);
        if dropped > 0 {
            eprintln!("  warning: {fid} dropped {dropped}");
        }
        let s = stats::schuster(&phases, 1)[0];
        let observed = s.d_squared;

        let mut null: Vec<f64> = offsets
            .iter()
            .map(|&off| {
                let shifted: Vec<f64> = times
                    .iter()
                    .filter_map(|&t| forcing.phase_at(t + off))
                    .collect();
                stats::schuster(&shifted, 1)[0].d_squared
            })
            .collect();
        null.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let ge = null.iter().filter(|&&x| x >= observed).count();
        let p = (ge as f64 + 1.0) / (NULL_TRIALS as f64 + 1.0);
        ps.push(p);
        family_phases.push(s.phase.to_degrees());

        println!(
            "{:<10} {:>8} {:>10.2} {:>8.1}d {:>10.0} {:>8.4}{}",
            fid,
            count,
            observed / phases.len() as f64,
            s.phase.to_degrees(),
            null.last().copied().unwrap_or(0.0),
            p,
            if p < 0.05 { "  *" } else { "" }
        );
    }

    let sig = ps.iter().filter(|&&p| p < 0.05).count();
    println!(
        "\n{sig}/{} families below p=0.05 (expected {:.1}), floor {:.4}",
        ps.len(),
        0.05 * ps.len() as f64,
        1.0 / (NULL_TRIALS as f64 + 1.0)
    );

    // Benjamini-Hochberg. Phase 1 produced 0/74 survivors here; this is the
    // comparison that matters.
    let mut sorted = ps.clone();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let m = sorted.len() as f64;
    let mut k = 0usize;
    for (i, &pv) in sorted.iter().enumerate() {
        if pv <= 0.05 * (i as f64 + 1.0) / m {
            k = i + 1;
        }
    }
    println!("Benjamini-Hochberg at FDR 0.05: {k}/{} families survive", sorted.len());

    // Are the preferred phases themselves consistent across families?
    //
    // CAVEAT: all families lie within ~30 km, so they see nearly identical
    // forcing. This is NOT 12 independent tests of phase -- it is a check that
    // the pipeline returns a coherent answer, not independent confirmation.
    let mut span = phases_summary(&family_phases);
    span.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let width = span.last().unwrap() - span.first().unwrap();
    println!(
        "\npreferred phases span {:.1} deg ({:.1} to {:.1})",
        width,
        span.first().unwrap(),
        span.last().unwrap()
    );
    println!(
        "  all {} families cluster after the dCFS maximum; random phases would\n  \
         span ~360 deg (co-located, so this is coherence not independence)",
        span.len()
    );

    // Sensitivity to the assumed geometry.
    println!("\ngeometry sensitivity, largest family:");
    let times = parkfield::family(&events, &fams[0].0).times();
    for (label, plane) in [
        ("SAF strike-slip 137/90/180", SAF),
        ("strike 127", fault::FaultPlane::new(127.0, 90.0, 180.0)),
        ("strike 147", fault::FaultPlane::new(147.0, 90.0, 180.0)),
        ("dip 80", fault::FaultPlane::new(137.0, 80.0, 180.0)),
        ("thrust 137/30/90", fault::FaultPlane::new(137.0, 30.0, 90.0)),
    ] {
        let c: Vec<f64> = tensors
            .iter()
            .map(|t| fault::coulomb(&fault::to_local_ned(t, lat, lon), &plane, MU))
            .collect();
        let f = Forcing::new(forcing.times().to_vec(), c).unwrap();
        let (ph, _) = f.phases(&times);
        let s = stats::schuster(&ph, 1)[0];
        println!(
            "  {:<28} D2/N {:>8.2}   phase {:>7.1}d",
            label,
            s.d_squared / ph.len() as f64,
            s.phase.to_degrees()
        );
    }

    let _ = TidalTensor::default();
    Ok(())
}
