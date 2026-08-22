//! A1: recompute the phase lags with **total** tide — solid plus ocean loading.
//!
//!     cargo run --release --example total_tide_lag
//!
//! # The prediction, fixed before running
//!
//! H1 found lags consistent *within* each band but 132.6° apart *between* them:
//! semidiurnal −9.4° ± 12.1°, diurnal +123.2° ± 1.6°. Two readings were offered:
//! a real frequency-dependent phase response, or a **band-dependent error in the
//! ΔCFS reference** caused by omitting ocean loading, whose phase relative to the
//! solid tide differs between bands.
//!
//! > **If loading is the cause, adding it brings the two bands' lags together.**
//! > **If they stay ~130° apart, the frequency-dependent response is real.**
//!
//! Either outcome is informative. This also retires the circular reasoning noted
//! in the self-critique: the flat R(ω) that killed the band prediction's premise
//! was itself computed with the forcing term this run restores.
//!
//! # Loading constants
//!
//! From SPOTL/GOT4.7 at each site's centroid, in nanostrain, as
//! `(amplitude, phase)` for east extension, north extension and tensor shear.
//! Reconstructed as
//! `ε(t) = A cos(χ_local(t) + φ + 10.97°)`, the offset verified for M2 in
//! `verify_loading_phase.rs`.
//!
//! ⚠ **The 10.97° calibration is verified for M2 only** — A4 showed `hartid`'s
//! constituent inference blocks cross-band validation. Applying it to N2, O1 and
//! Q1 is an assumption, and a wrong per-constituent offset would move the very
//! lags being measured. Stated here rather than buried.
//!
//! ⚠ Cascadia tremor spans 40–50 °N; one centroid is a crude stand-in for a field
//! that varies sharply along a coastline.

use eqf::{cascadia, parkfield};
use ph_core::{doodson, fault, field::TidalField, love::Elastic};
use rustspice_core::{Et, KernelSet};

const PHASE_OFFSET_DEG: f64 = 10.97;
const STEP: f64 = 0.02;
const BANDS: &[&str] = &["M2", "N2", "O1", "Q1"];

/// `(ee_amp, ee_pha, nn_amp, nn_pha, en_amp, en_pha)` in nanostrain and degrees.
type Load = [f64; 6];

struct Site {
    name: &'static str,
    lat: f64,
    lon: f64,
    plane: fault::FaultPlane,
    /// Loading per constituent, in the order of `BANDS`.
    load: [Load; 4],
}

fn sites() -> Vec<Site> {
    vec![
        Site {
            name: "Parkfield",
            lat: 35.635,
            lon: -120.150,
            // Deep San Andreas: vertical right-lateral strike-slip.
            plane: fault::FaultPlane::new(137.0, 90.0, 180.0),
            load: [
                [3.8812, 25.2010, 3.4624, 145.3707, 2.7008, 130.4305],
                [0.7939, 47.9213, 0.6804, 165.1476, 0.6156, 150.1291],
                [0.7289, -110.1990, 0.4834, 19.9494, 0.8818, -62.4000],
                [0.1241, -110.5112, 0.0888, 15.9653, 0.1483, -50.8124],
            ],
        },
        Site {
            name: "Cascadia",
            lat: 44.541,
            lon: -123.422,
            // Subduction megathrust: shallow-dipping thrust.
            plane: fault::FaultPlane::new(350.0, 12.0, 90.0),
            load: [
                [11.2122, 15.7136, 6.8562, -175.1720, 4.0784, 166.4862],
                [2.3373, 39.6965, 1.3852, -150.2765, 0.8707, -172.0249],
                [2.8839, -104.4336, 2.1259, 74.4423, 0.6735, 32.8543],
                [0.5106, -100.7238, 0.3717, 79.2016, 0.1350, 37.0955],
            ],
        },
    ]
}

/// Loading stress tensor in local NED, from horizontal strains under a free
/// surface (`σ_zz = 0`, plane stress, μ = 30 GPa, ν = 0.25).
fn loading_tensor(site: &Site, day: f64) -> ph_core::tidal::TidalTensor {
    const NU: f64 = 0.25;
    const MU: f64 = 3.0e10;
    let e_mod = 2.0 * MU * (1.0 + NU);
    let off = PHASE_OFFSET_DEG.to_radians();

    let (mut ee, mut nn, mut en) = (0.0f64, 0.0f64, 0.0f64);
    for (i, name) in BANDS.iter().enumerate() {
        let c = doodson::constituent(name).unwrap();
        let chi = c.phase_at_longitude(day, site.lon);
        let l = &site.load[i];
        // nanostrain -> strain
        ee += l[0] * 1e-9 * (chi + l[1].to_radians() + off).cos();
        nn += l[2] * 1e-9 * (chi + l[3].to_radians() + off).cos();
        en += l[4] * 1e-9 * (chi + l[5].to_radians() + off).cos();
    }

    let s_ee = e_mod / (1.0 - NU * NU) * (ee + NU * nn);
    let s_nn = e_mod / (1.0 - NU * NU) * (nn + NU * ee);
    let s_en = e_mod / (1.0 + NU) * en;
    // NED order; the free surface gives sigma_DD = 0.
    ph_core::tidal::TidalTensor {
        m: [[s_nn, s_en, 0.0], [s_en, s_ee, 0.0], [0.0, 0.0, 0.0]],
    }
}

fn norm180(d: f64) -> f64 {
    let mut x = d % 360.0;
    if x > 180.0 {
        x -= 360.0;
    }
    if x < -180.0 {
        x += 360.0;
    }
    x
}

fn main() -> rustspice_core::Result<()> {
    let mut ks = KernelSet::new();
    for k in ["naif0012.tls", "de440s.bsp", "pck00011.tpc", "gm_de440.tpc"] {
        ks.add_file(format!("kernels/{k}"))?;
    }
    let mut spice = ks.open()?;
    let elastic = Elastic::EARTH;
    let field = TidalField::on_earth(&mut spice, "IAU_EARTH")?;

    let pk: Vec<f64> = {
        let e = parkfield::parse_catalog(
            &std::fs::read_to_string("data/parkfield/LFEcat_Apr2001-Apr2024.csv").expect("pk"),
        );
        let mut t: Vec<f64> = e.iter().map(|x| x.day).collect();
        t.sort_by(|a, b| a.partial_cmp(b).unwrap());
        t
    };
    let cs: Vec<f64> = {
        let e = cascadia::parse_catalog(
            &std::fs::read_to_string("data/cascadia/cascadia_tremor.csv").expect("cs"),
        );
        let mut t: Vec<f64> = e.iter().map(|x| x.day).collect();
        t.sort_by(|a, b| a.partial_cmp(b).unwrap());
        t
    };

    for (site, times) in sites().into_iter().zip([pk, cs]) {
        let (t0, t1) = (times[0], times[times.len() - 1]);
        let n = ((t1 - t0) / STEP) as usize;
        let epoch2000 = spice.parse_time("2000-01-01T00:00:00")?;
        let days: Vec<f64> = (0..n).map(|i| t0 + i as f64 * STEP).collect();
        let epochs: Vec<Et> = days.iter().map(|&d| Et(epoch2000.0 + d * 86400.0)).collect();
        let tensors = field.tensors(&mut spice, &epochs)?;

        // Solid-only and total dCFS series.
        let solid: Vec<f64> = tensors
            .iter()
            .map(|t| {
                elastic.stress(fault::coulomb(
                    &fault::to_local_ned(t, site.lat, site.lon),
                    &site.plane,
                    0.4,
                ))
            })
            .collect();
        let total: Vec<f64> = tensors
            .iter()
            .zip(&days)
            .map(|(t, &d)| {
                let mut m = fault::to_local_ned(t, site.lat, site.lon);
                let s = elastic.stress_per_tensor();
                let load = loading_tensor(&site, d);
                for i in 0..3 {
                    for j in 0..3 {
                        // Solid tensor is in s^-2; convert, then add loading in Pa.
                        m.m[i][j] = m.m[i][j] * s + load.m[i][j] / s * s;
                    }
                }
                // m is now in Pa throughout.
                fault::coulomb(&m, &site.plane, 0.4)
            })
            .collect();

        println!("\n=== {} : {} events ===", site.name, times.len());
        println!(
            "{:<5} {:>10} {:>10} {:>10} {:>10} {:>9}",
            "band", "solid_amp", "total_amp", "lag_solid", "lag_total", "change"
        );
        println!("{}", "-".repeat(58));

        let mut solid_lags = Vec::new();
        let mut total_lags = Vec::new();
        for name in BANDS {
            let c = doodson::constituent(name).unwrap();
            let fit = |series: &[f64]| -> (f64, f64) {
                let (mut a, mut b) = (0.0f64, 0.0f64);
                for (k, &d) in days.iter().enumerate() {
                    let (s, co) = c.phase_at_longitude(d, site.lon).sin_cos();
                    a += series[k] * co;
                    b += series[k] * s;
                }
                (
                    2.0 * (a * a + b * b).sqrt() / days.len() as f64,
                    (-b).atan2(a).to_degrees(),
                )
            };
            let (amp_s, pha_s) = fit(&solid);
            let (amp_t, pha_t) = fit(&total);

            let (mut ea, mut eb) = (0.0f64, 0.0f64);
            for &t in &times {
                let (s, co) = c.phase_at_longitude(t, site.lon).sin_cos();
                ea += co;
                eb += s;
            }
            let obs = eb.atan2(ea).to_degrees();
            let ls = norm180(obs - pha_s);
            let lt = norm180(obs - pha_t);
            solid_lags.push((c.period_days(), ls));
            total_lags.push((c.period_days(), lt));
            println!(
                "{:<5} {:>10.1} {:>10.1} {:>10.1} {:>10.1} {:>9.1}",
                name,
                amp_s,
                amp_t,
                ls,
                lt,
                norm180(lt - ls)
            );
        }

        let band = |v: &[(f64, f64)], semi: bool| -> f64 {
            let s: Vec<f64> = v
                .iter()
                .filter(|(p, _)| (*p < 0.75) == semi)
                .map(|(_, l)| *l)
                .collect();
            s.iter().sum::<f64>() / s.len() as f64
        };
        let ds = norm180(band(&solid_lags, false) - band(&solid_lags, true));
        let dt = norm180(band(&total_lags, false) - band(&total_lags, true));
        println!(
            "  band difference (diurnal - semidiurnal):  solid {:+.1} deg  ->  total {:+.1} deg",
            ds, dt
        );
        println!(
            "  {}",
            if dt.abs() < ds.abs() * 0.5 {
                "CONVERGED -- loading explains the band split; the prediction holds"
            } else if dt.abs() < ds.abs() * 0.9 {
                "partial convergence -- loading explains some of the split"
            } else {
                "NO convergence -- the band split is not a loading artifact"
            }
        );
    }
    Ok(())
}
