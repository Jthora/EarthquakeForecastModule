//! C4 (second attempt): response per constituent, from analytic Doodson phases.
//!
//!     cargo run --release --example parkfield_constituent_response
//!
//! # What changed
//!
//! The first attempt derived constituent phase by demodulating the composite
//! ΔCFS. That failed twice (traps 5 and 6). Phase now comes from
//! [`ph_core::doodson`], where each constituent's argument is an integer
//! combination of the six fundamental astronomical arguments. Its phase is
//! **uniform over long spans by construction** — verified in that module's tests
//! to better than 5% across 12 histogram bins.
//!
//! # The null: per-block random shifts
//!
//! A **global** time shift cannot work here and never will. `D²` is invariant
//! under rotation, and shifting a single constituent's phase globally *is* a
//! rotation — concentration is untouched. That is trap 5, and it is a property of
//! the statistic, not a fixable detail.
//!
//! So the null shifts **each block independently**. Blocks are long compared with
//! the constituent period, so within-block clustering (aftershock-like bursts,
//! detection outages) is preserved, while the alignment between blocks is
//! randomised. That changes concentration, which is what a null must do.
//!
//! Block length is `max(4 × period, 30 d)`. Long-period constituents therefore get
//! few blocks, and the block count is reported so an underpowered test is visible
//! rather than silent.

use eqf::{parkfield};
use ph_core::{doodson, fault, field::TidalField, love::Elastic};
use rustspice_core::{Et, KernelSet};

const SAF: fault::FaultPlane = fault::FaultPlane {
    strike_deg: 137.0,
    dip_deg: 90.0,
    rake_deg: 180.0,
};

const NULL_TRIALS: usize = 400;
const USABLE: &[&str] = &["M2", "N2", "O1", "Q1", "Mf", "Msf", "Mm", "Ssa", "Sa"];

struct Rng(u64);
impl Rng {
    fn next_f64(&mut self) -> f64 {
        self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        ((self.0 >> 11) as f64) / ((1u64 << 53) as f64)
    }
}

fn schuster_power(phases: &[f64]) -> f64 {
    let (mut a, mut b) = (0.0f64, 0.0f64);
    for &p in phases {
        let (s, c) = p.sin_cos();
        a += c;
        b += s;
    }
    (a * a + b * b) / phases.len() as f64
}

fn main() -> rustspice_core::Result<()> {
    let events = parkfield::parse_catalog(
        &std::fs::read_to_string("data/parkfield/LFEcat_Apr2001-Apr2024.csv")
            .expect("run scripts/fetch-parkfield.sh"),
    );
    let mut times: Vec<f64> = events.iter().map(|e| e.day).collect();
    times.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let (t0, t1) = (times[0], times[times.len() - 1]);
    println!("{} events, span {:.0} d", times.len(), t1 - t0);

    // Constituent AMPLITUDES in our own dCFS, by harmonic regression on the
    // analytic Doodson phase. This is the right normalisation: not the raw tidal
    // potential, but the Coulomb stress this fault actually sees at this site.
    //
    // Without it a falling D2/N with period cannot be read as a transfer
    // function, because the tidal potential's own amplitudes fall with period
    // too. Response per unit stress is the quantity that separates them.
    let mut ks = KernelSet::new();
    for k in ["naif0012.tls", "de440s.bsp", "pck00011.tpc", "gm_de440.tpc"] {
        ks.add_file(format!("kernels/{k}"))?;
    }
    let mut spice = ks.open()?;
    let (lat, lon, _) = parkfield::family_location(
        &events,
        &parkfield::families(&events)[0].0,
    )
    .unwrap();
    let step = 0.02;
    let n = ((t1 - t0) / step) as usize;
    let epoch2000 = spice.parse_time("2000-01-01T00:00:00")?;
    let sdays: Vec<f64> = (0..n).map(|i| t0 + i as f64 * step).collect();
    let epochs: Vec<Et> = sdays.iter().map(|&d| Et(epoch2000.0 + d * 86400.0)).collect();
    let earth = TidalField::on_earth(&mut spice, "IAU_EARTH")?;
    let tensors = earth.tensors(&mut spice, &epochs)?;
    let elastic = Elastic::EARTH;
    let cfs: Vec<f64> = tensors
        .iter()
        .map(|t| elastic.stress(fault::coulomb(&fault::to_local_ned(t, lat, lon), &SAF, 0.4)))
        .collect();
    println!("dCFS sampled at {n} points, converted to Pa via Love numbers\n");

    println!(
        "{:<6} {:>9} {:>9} {:>10} {:>9} {:>10} {:>8}",
        "band", "period", "amp(Pa)", "D2/N", "ratio", "R(w)", "p"
    );
    println!("{}", "-".repeat(70));

    let mut rng = Rng(0x0D00D5);
    let mut results = Vec::new();

    for name in USABLE {
        let c = doodson::constituent(name).unwrap();
        let period = c.period_days();
        let block = (4.0 * period).max(30.0);
        let n_blocks = ((t1 - t0) / block).floor() as usize;

        let observed = schuster_power(&c.phases(&times));

        // Least-squares amplitude of this constituent in the dCFS series.
        let (mut sa, mut sb) = (0.0f64, 0.0f64);
        for (k, &d) in sdays.iter().enumerate() {
            let (si, co) = c.phase_at(d).sin_cos();
            sa += cfs[k] * co;
            sb += cfs[k] * si;
        }
        let amp = 2.0 * (sa * sa + sb * sb).sqrt() / sdays.len() as f64;

        let mut null: Vec<f64> = (0..NULL_TRIALS)
            .map(|_| {
                // One independent offset per block, drawn from a full cycle.
                let offs: Vec<f64> = (0..=n_blocks + 1)
                    .map(|_| rng.next_f64() * period)
                    .collect();
                let shifted: Vec<f64> = times
                    .iter()
                    .map(|&t| {
                        let b = (((t - t0) / block).floor() as usize).min(offs.len() - 1);
                        c.phase_at(t + offs[b])
                    })
                    .collect();
                schuster_power(&shifted)
            })
            .collect();
        null.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let ge = null.iter().filter(|&&x| x >= observed).count();
        let p = (ge as f64 + 1.0) / (NULL_TRIALS as f64 + 1.0);

        // The null median rises with period, because longer blocks mean fewer
        // independent randomisations. Raw D2/N is therefore NOT comparable across
        // constituents; the ratio to each band's own null median is.
        let med = null[null.len() / 2];
        // Fractional rate modulation, then response per unit stress.
        let eps = 2.0 * (observed / times.len() as f64).sqrt();
        println!(
            "{:<6} {:>9.4} {:>9.2} {:>10.1} {:>9.1} {:>10.3e} {:>8.4}{}{}",
            name,
            period,
            amp,
            observed,
            observed / med,
            eps / amp,
            p,
            if p < 0.05 { "  *" } else { "" },
            if n_blocks < 8 { "  [few]" } else { "" }
        );
        results.push((name.to_string(), p, n_blocks));
    }

    // Benjamini-Hochberg over the constituents tested.
    let mut ps: Vec<f64> = results.iter().map(|(_, p, _)| *p).collect();
    ps.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let m = ps.len() as f64;
    let mut k = 0usize;
    for (i, &pv) in ps.iter().enumerate() {
        if pv <= 0.05 * (i as f64 + 1.0) / m {
            k = i + 1;
        }
    }
    println!("\nBenjamini-Hochberg at FDR 0.05: {k}/{} constituents survive", ps.len());
    println!("null floor 1/(n+1) = {:.4}", 1.0 / (NULL_TRIALS as f64 + 1.0));
    println!("\namp  = dCFS amplitude at this constituent, Pa (Love-calibrated, ~2x)");
    println!("R(w) = eps/amp: fractional rate modulation per Pa -- the transfer function");
    println!("\nK1, S1, S2, P1, K2 excluded: locked to solar time and degenerate");
    println!("with the diurnal/thermal detection artifact measured in D1.");
    Ok(())
}
