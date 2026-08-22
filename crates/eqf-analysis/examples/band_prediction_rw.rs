//! The band prediction test: R(ω) with total tide, from M2 to Ssa.
//!
//!     cargo run --release --example band_prediction_rw
//!
//! # What makes this the real test
//!
//! Every earlier R(ω) was solid-tide only, and the flat result was used to kill the
//! band prediction's premise — circularly, since the missing term was deferred on
//! the strength of that same flat result.
//!
//! Ocean loading is now included at all seven constituents, and it is strongly
//! **band-dependent**: at Parkfield it changes M2 forcing by −28% and O1 by +26%,
//! while Mf loading is 1.5% of M2's and Ssa is essentially zero. So the correction
//! reshapes the short-period end and leaves the long end intact — exactly the
//! asymmetry needed to test whether R rises toward long periods.
//!
//! # Bounds, not point estimates
//!
//! Long-period constituents are **not significant** at either site. R = ε/forcing
//! is meaningless there as a point estimate — ε is noise divided by a small number,
//! the error that made Sa's R come out 20× everything else earlier.
//!
//! So the long end is reported as an **upper bound** from the null maximum:
//!
//! | Outcome | Verdict |
//! |---|---|
//! | `R_upper(long) < R(M2)` | band prediction **refuted** |
//! | `R_upper(long) ≫ R(M2)` | **unconstrained** — say so |
//!
//! # Geometry
//!
//! Cascadia's fault plane was invented in the previous run. Rather than substitute
//! false precision, three plausible megathrust orientations are run and the
//! **shape** of R(ω) compared. If the shape is robust to geometry, the assumption
//! does not matter for this question.

use eqf::{cascadia, parkfield};
use ph_core::{doodson, fault, field::TidalField, love::Elastic, tidal::TidalTensor};
use rustspice_core::{Et, KernelSet};

const PHASE_OFFSET_DEG: f64 = 10.97;
const STEP: f64 = 0.05;
const NULL_TRIALS: usize = 200;
const BANDS: &[&str] = &["M2", "N2", "O1", "Q1", "Mf", "Mm", "Ssa"];

/// (ee_amp, ee_pha, nn_amp, nn_pha, en_amp, en_pha) nanostrain, degrees.
const PK_LOAD: [[f64; 6]; 7] = [
    [3.8956, 24.39, 3.4984, 145.17, 2.6620, 132.55],
    [0.7584, 46.50, 0.6699, 164.19, 0.5956, 153.34],
    [0.7448, -110.22, 0.4954, 19.46, 0.8549, -61.03],
    [0.1258, -108.29, 0.0909, 19.37, 0.1471, -50.80],
    [0.0588, -160.70, 0.0385, -2.51, 0.0319, -33.02],
    [0.0218, -160.61, 0.0198, 6.04, 0.0196, -23.20],
    [0.0169, -180.00, 0.0198, 0.00, 0.0208, 0.00],
];
const CAS_LOAD: [[f64; 6]; 7] = [
    [11.0051, 15.81, 6.6519, -175.20, 4.0509, 167.90],
    [2.2683, 39.93, 1.3237, -149.89, 0.8633, -169.01],
    [2.8530, -104.93, 2.0653, 74.14, 0.7179, 33.44],
    [0.4993, -99.42, 0.3610, 80.07, 0.1338, 41.86],
    [0.1674, -158.22, 0.1422, 12.34, 0.0737, 3.20],
    [0.0628, -164.43, 0.0634, 8.04, 0.0377, -2.70],
    [0.0440, -180.00, 0.0506, 0.00, 0.0314, 0.00],
];

struct Rng(u64);
impl Rng {
    fn next_f64(&mut self) -> f64 {
        self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        ((self.0 >> 11) as f64) / ((1u64 << 53) as f64)
    }
}

fn loading_tensor(load: &[[f64; 6]; 7], lon: f64, day: f64) -> TidalTensor {
    const NU: f64 = 0.25;
    const MU: f64 = 3.0e10;
    let e_mod = 2.0 * MU * (1.0 + NU);
    let off = PHASE_OFFSET_DEG.to_radians();
    let (mut ee, mut nn, mut en) = (0.0f64, 0.0f64, 0.0f64);
    for (i, name) in BANDS.iter().enumerate() {
        let chi = doodson::constituent(name).unwrap().phase_at_longitude(day, lon);
        let l = &load[i];
        ee += l[0] * 1e-9 * (chi + l[1].to_radians() + off).cos();
        nn += l[2] * 1e-9 * (chi + l[3].to_radians() + off).cos();
        en += l[4] * 1e-9 * (chi + l[5].to_radians() + off).cos();
    }
    let s_ee = e_mod / (1.0 - NU * NU) * (ee + NU * nn);
    let s_nn = e_mod / (1.0 - NU * NU) * (nn + NU * ee);
    let s_en = e_mod / (1.0 + NU) * en;
    TidalTensor { m: [[s_nn, s_en, 0.0], [s_en, s_ee, 0.0], [0.0, 0.0, 0.0]] }
}

#[allow(clippy::too_many_arguments)]
fn run(
    label: &str,
    times: &[f64],
    lat: f64,
    lon: f64,
    plane: fault::FaultPlane,
    load: &[[f64; 6]; 7],
    tensors: &[TidalTensor],
    days: &[f64],
    rng: &mut Rng,
    quiet: bool,
) -> Vec<(f64, f64, bool)> {
    let e = Elastic::EARTH;
    let scale = e.stress_per_tensor();
    let cfs: Vec<f64> = tensors
        .iter()
        .zip(days)
        .map(|(t, &d)| {
            let mut m = fault::to_local_ned(t, lat, lon);
            let l = loading_tensor(load, lon, d);
            for i in 0..3 {
                for j in 0..3 {
                    m.m[i][j] = m.m[i][j] * scale + l.m[i][j];
                }
            }
            fault::coulomb(&m, &plane, 0.4)
        })
        .collect();

    if !quiet {
        println!("\n=== {label} ===");
        println!(
            "{:<5} {:>9} {:>11} {:>10} {:>11} {:>7}",
            "band", "period_d", "force_Pa", "eps", "R or bound", "p"
        );
        println!("{}", "-".repeat(58));
    }

    let mut out = Vec::new();
    for (_i, name) in BANDS.iter().enumerate() {
        let c = doodson::constituent(name).unwrap();
        let (mut a, mut b) = (0.0f64, 0.0f64);
        for (k, &d) in days.iter().enumerate() {
            let (s, co) = c.phase_at_longitude(d, lon).sin_cos();
            a += cfs[k] * co;
            b += cfs[k] * s;
        }
        let force = 2.0 * (a * a + b * b).sqrt() / days.len() as f64;

        let pw = |shift: &dyn Fn(f64) -> f64| {
            let (mut x, mut y) = (0.0f64, 0.0f64);
            for &t in times {
                let (s, co) = c.phase_at_longitude(t + shift(t), lon).sin_cos();
                x += co;
                y += s;
            }
            (x * x + y * y) / times.len() as f64
        };
        let obs = pw(&|_| 0.0);
        let period = c.period_days();
        let block = (4.0 * period).max(30.0);
        let (t0, t1) = (times[0], times[times.len() - 1]);
        let nb = (((t1 - t0) / block).floor() as usize) + 2;
        let mut null: Vec<f64> = (0..NULL_TRIALS)
            .map(|_| {
                let offs: Vec<f64> = (0..nb).map(|_| rng.next_f64() * period).collect();
                pw(&|t| offs[(((t - t0) / block).floor() as usize).min(offs.len() - 1)])
            })
            .collect();
        null.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let ge = null.iter().filter(|&&x| x >= obs).count();
        let p = (ge as f64 + 1.0) / (NULL_TRIALS as f64 + 1.0);
        let sig = p < 0.05;

        let eps = 2.0 * (obs / times.len() as f64).sqrt();
        let eps_bound = 2.0 * (null[null.len() - 1] / times.len() as f64).sqrt();
        let r = if sig { eps / force } else { eps_bound / force };
        if !quiet {
            println!(
                "{:<5} {:>9.3} {:>11.2} {:>9.2}% {:>10.2e}{} {:>7.4}",
                name,
                period,
                force,
                eps * 100.0,
                r,
                if sig { " " } else { "<" },
                p
            );
        }
        out.push((period, r, sig));
    }
    out
}

fn verdict(rows: &[(f64, f64, bool)]) {
    let m2 = rows[0].1;
    let long: Vec<&(f64, f64, bool)> = rows.iter().filter(|r| r.0 > 5.0).collect();
    let worst = long.iter().map(|r| r.1).fold(0.0, f64::max);
    println!(
        "  R(M2) = {m2:.2e};  long-period bound reaches {worst:.2e}  ->  ratio {:.1}",
        worst / m2
    );
    println!(
        "  {}",
        if worst < m2 {
            "long-period response is BOUNDED BELOW R(M2) -- band prediction REFUTED"
        } else {
            "long-period bound exceeds R(M2) -- UNCONSTRAINED, the prediction survives untested"
        }
    );
}

fn main() -> rustspice_core::Result<()> {
    let mut ks = KernelSet::new();
    for k in ["naif0012.tls", "de440s.bsp", "pck00011.tpc", "gm_de440.tpc"] {
        ks.add_file(format!("kernels/{k}"))?;
    }
    let mut spice = ks.open()?;
    let field = TidalField::on_earth(&mut spice, "IAU_EARTH")?;
    let mut rng = Rng(0xBA4D2);

    for (label, path, lat, lon, load, planes) in [
        (
            "Parkfield",
            "data/parkfield/LFEcat_Apr2001-Apr2024.csv",
            35.635,
            -120.150,
            &PK_LOAD,
            vec![("SAF 137/90/180", fault::FaultPlane::new(137.0, 90.0, 180.0))],
        ),
        (
            "Cascadia",
            "data/cascadia/cascadia_tremor.csv",
            44.541,
            -123.422,
            &CAS_LOAD,
            vec![
                ("megathrust 350/12/90", fault::FaultPlane::new(350.0, 12.0, 90.0)),
                ("megathrust 353/18/90", fault::FaultPlane::new(353.0, 18.0, 90.0)),
                ("megathrust 000/25/90", fault::FaultPlane::new(0.0, 25.0, 90.0)),
            ],
        ),
    ] {
        let txt = std::fs::read_to_string(path).expect("catalogue");
        let mut times: Vec<f64> = if label == "Parkfield" {
            parkfield::parse_catalog(&txt).iter().map(|e| e.day).collect()
        } else {
            cascadia::parse_catalog(&txt).iter().map(|e| e.day).collect()
        };
        times.sort_by(|a, b| a.partial_cmp(b).unwrap());

        let (t0, t1) = (times[0], times[times.len() - 1]);
        let n = ((t1 - t0) / STEP) as usize;
        let epoch2000 = spice.parse_time("2000-01-01T00:00:00")?;
        let days: Vec<f64> = (0..n).map(|i| t0 + i as f64 * STEP).collect();
        let epochs: Vec<Et> = days.iter().map(|&d| Et(epoch2000.0 + d * 86400.0)).collect();
        let tensors = field.tensors(&mut spice, &epochs)?;

        for (i, (gname, plane)) in planes.iter().enumerate() {
            let rows = run(
                &format!("{label} : {gname}"),
                &times,
                lat,
                lon,
                *plane,
                load,
                &tensors,
                &days,
                &mut rng,
                false,
            );
            verdict(&rows);
            if i == 0 && planes.len() > 1 {
                println!("  (geometry sensitivity follows)");
            }
        }
    }
    Ok(())
}
