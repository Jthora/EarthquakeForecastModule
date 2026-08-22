//! P3.2: does the tide promote failure on the faults that actually broke?
//!
//!     cargo run --release --example gcmt_coulomb_sign
//!
//! # Why this is a different test, not another slice
//!
//! P3.4 and P3.5 used raw tidal phase and found nothing, bounding ordinary crust
//! below ~4%. But raw phase is **blind to whether the tide loads or unloads each
//! fault**: a tide promoting failure on a thrust in Japan simultaneously inhibits
//! it on a normal fault in Greece, and pooled globally those cancel.
//!
//! With GCMT mechanisms we can resolve ΔCFS onto each event's *own* fault plane
//! and ask a question raw phase cannot:
//!
//! > **Do earthquakes occur preferentially when ΔCFS on their own fault is
//! > positive — that is, when the tide is promoting failure?**
//!
//! Under no effect the answer is 50%. This raises signal by aligning the feature
//! with the physics rather than by discarding events, which is why it is not the
//! stratification P3.5 showed makes bounds worse.
//!
//! # Nodal plane ambiguity
//!
//! A moment tensor gives two equally valid planes and does not say which broke.
//! Both are run separately as a pre-registered robustness check. Using the wrong
//! plane for roughly half the events dilutes a real signal; it cannot manufacture
//! one.
//!
//! # Null
//!
//! Block-shift, as elsewhere: aftershock sequences make event times strongly
//! non-independent, so the binomial 50% is not a valid null on its own.

use eqf::gcmt;
use ph_core::{fault, field::TidalField, love::Elastic};
use rustspice_core::{Et, KernelSet};

const MIN_MW: f64 = 5.5;
const MU: f64 = 0.4;
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

    let all = gcmt::parse_ndk(
        &std::fs::read_to_string("data/gcmt/gcmt.ndk").expect("run scripts/fetch-gcmt.sh"),
    );
    let ev = gcmt::select(&all, MIN_MW, None);
    println!(
        "GCMT: {} solutions, {} at Mw >= {MIN_MW}, {:.1} yr span",
        all.len(),
        ev.len(),
        (ev[ev.len() - 1].day - ev[0].day) / 365.25
    );

    let epoch2000 = spice.parse_time("2000-01-01T00:00:00")?;
    let earth = TidalField::on_earth(&mut spice, "IAU_EARTH")?;
    let elastic = Elastic::EARTH;

    // One tensor evaluation per event, at the event's own centroid time.
    let epochs: Vec<Et> = ev.iter().map(|e| Et(epoch2000.0 + e.day * 86400.0)).collect();
    println!("computing tidal tensors at {} event times...", epochs.len());
    let tensors = earth.tensors(&mut spice, &epochs)?;

    // Local NED tensors, reused across both planes and every null trial.
    let local: Vec<_> = tensors
        .iter()
        .zip(&ev)
        .map(|(t, e)| fault::to_local_ned(t, e.lat_deg, e.lon_deg))
        .collect();

    // For the null we need tensors at shifted times too. Precompute a shared
    // shifted set per trial would be huge, so instead evaluate the null by
    // shifting into a precomputed dense grid is impractical globally -- events
    // are at different longitudes. Evaluate directly per trial instead.
    let t0 = ev[0].day;

    let mut rng = Rng(0x6C47);
    for (which, label) in [(0usize, "nodal plane 1"), (1, "nodal plane 2")] {
        let cfs: Vec<f64> = local
            .iter()
            .zip(&ev)
            .map(|(t, e)| elastic.stress(fault::coulomb(t, &e.planes()[which], MU)))
            .collect();

        let positive = cfs.iter().filter(|&&c| c > 0.0).count();
        let frac = positive as f64 / cfs.len() as f64;
        let mean = cfs.iter().sum::<f64>() / cfs.len() as f64;

        // Null: shift each block by a random offset drawn from a lunar day, then
        // recompute. Requires fresh tensors, so this is the expensive part.
        let n_blocks = (((ev[ev.len() - 1].day - t0) / BLOCK_DAYS).floor() as usize) + 2;
        let mut null_frac = Vec::with_capacity(NULL_TRIALS);
        let mut null_mean = Vec::with_capacity(NULL_TRIALS);
        for _ in 0..NULL_TRIALS {
            let offs: Vec<f64> = (0..n_blocks).map(|_| rng.next_f64() * 1.03505).collect();
            let shifted: Vec<Et> = ev
                .iter()
                .map(|e| {
                    let b = (((e.day - t0) / BLOCK_DAYS).floor() as usize).min(offs.len() - 1);
                    Et(epoch2000.0 + (e.day + offs[b]) * 86400.0)
                })
                .collect();
            let ts = earth.tensors(&mut spice, &shifted)?;
            let c: Vec<f64> = ts
                .iter()
                .zip(&ev)
                .map(|(t, e)| {
                    let l = fault::to_local_ned(t, e.lat_deg, e.lon_deg);
                    elastic.stress(fault::coulomb(&l, &e.planes()[which], MU))
                })
                .collect();
            null_frac.push(c.iter().filter(|&&x| x > 0.0).count() as f64 / c.len() as f64);
            null_mean.push(c.iter().sum::<f64>() / c.len() as f64);
        }
        null_frac.sort_by(|a, b| a.partial_cmp(b).unwrap());
        null_mean.sort_by(|a, b| a.partial_cmp(b).unwrap());

        let ge_f = null_frac.iter().filter(|&&x| x >= frac).count();
        let ge_m = null_mean.iter().filter(|&&x| x >= mean).count();
        let p_f = (ge_f as f64 + 1.0) / (NULL_TRIALS as f64 + 1.0);
        let p_m = (ge_m as f64 + 1.0) / (NULL_TRIALS as f64 + 1.0);

        println!("\n=== {label} ===");
        println!(
            "  dCFS > 0 at:  {:.3}%   null median {:.3}%   p = {p_f:.4}{}",
            frac * 100.0,
            null_frac[null_frac.len() / 2] * 100.0,
            if p_f < 0.05 { "  *" } else { "" }
        );
        println!(
            "  mean dCFS  :  {mean:+.3} Pa   null median {:+.3} Pa   p = {p_m:.4}{}",
            null_mean[null_mean.len() / 2],
            if p_m < 0.05 { "  *" } else { "" }
        );
        println!(
            "  excess over 50%: {:+.3} percentage points",
            (frac - 0.5) * 100.0
        );
    }

    println!("\nnull floor 1/(n+1) = {:.4}", 1.0 / (NULL_TRIALS as f64 + 1.0));
    println!("A real effect should appear on at least one plane and not reverse sign");
    println!("between them; the true fault is one of the two for every event.");
    Ok(())
}
