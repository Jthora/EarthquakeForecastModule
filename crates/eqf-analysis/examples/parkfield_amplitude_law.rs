//! C3: does LFE response scale with tidal forcing amplitude?
//!
//!     cargo run --release --example parkfield_amplitude_law
//!
//! # Why this and not a constituent-by-constituent spectrum
//!
//! D1 showed the catalogue's diurnal and semidiurnal bands are dominated by a
//! detection artifact locked to solar time. Splitting response by constituent
//! runs straight back into that.
//!
//! But **the artifact does not care how strong the tide is.** So bin events by the
//! peak-to-trough amplitude of the ΔCFS cycle they fall in, and measure phase
//! concentration within each bin. Detection bias is amplitude-independent by
//! construction, so any monotonic trend is physical.
//!
//! This is a direct test of the amplitude law `R̃/r = Δτ/(aσ̄)` (Ader et al. 2014
//! eq. B7), and of whether Parkfield sits in the **non-linear** regime — Thomas et
//! al. (2012) infer `Aσ₀ = 6×10⁻⁴ MPa`, giving `S_T/Aσ₀ ≈ 0.2–2`, so
//! `R = R₀exp(S_T/Aσ₀)` should grow **faster than linearly** with amplitude.
//!
//! Predictions, fixed before the run:
//!
//! | Outcome | Reading |
//! |---|---|
//! | `D²/N` flat across amplitude bins | response is not tidal — artifact |
//! | `D²/N` grows ∝ amplitude² | linear response (`D` ∝ Δτ) |
//! | `D²/N` grows faster than amplitude² | non-linear, as the exponential predicts |

use eqf::{parkfield};
use ph_core::{fault, field::TidalField, phase::Forcing, stats};
use rustspice_core::{Et, KernelSet};

const SAF: fault::FaultPlane = fault::FaultPlane {
    strike_deg: 137.0,
    dip_deg: 90.0,
    rake_deg: 180.0,
};
const MU: f64 = 0.4;
const STEP_DAYS: f64 = 0.02;
const BINS: usize = 6;
const NULL_TRIALS: usize = 2000;

struct Rng(u64);
impl Rng {
    fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        self.0 >> 33
    }
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

    let (mut lo, mut hi) = (f64::MAX, f64::MIN);
    for e in &events {
        lo = lo.min(e.day);
        hi = hi.max(e.day);
    }
    let (t0, t1) = (lo - 4200.0, hi + 4200.0);
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
    let forcing = Forcing::new(days, cfs).expect("forcing");

    // Pool all families: the amplitude test needs counts, and every family sees
    // the same forcing anyway (they lie within ~30 km).
    let mut times: Vec<f64> = events.iter().map(|e| e.day).collect();
    times.sort_by(|a, b| a.partial_cmp(b).unwrap());

    // Pair each event with its cycle amplitude.
    let mut paired: Vec<(f64, f64)> = times
        .iter()
        .filter_map(|&t| forcing.cycle_amplitude_at(t).map(|a| (a, t)))
        .collect();
    paired.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    println!("{} events paired with a cycle amplitude\n", paired.len());

    let mut rng = Rng(0xA11B1AD);
    let offsets: Vec<f64> = (0..NULL_TRIALS)
        .map(|_| {
            let d = (30 + (rng.next_u64() % 3970) as i64) as f64;
            if rng.next_u64() % 2 == 0 { d } else { -d }
        })
        .collect();

    println!(
        "{:>4} {:>9} {:>12} {:>10} {:>10} {:>8}",
        "bin", "events", "amp(rel)", "D2/N", "null_max", "p"
    );
    println!("{}", "-".repeat(58));

    let per = paired.len() / BINS;
    let mut rows = Vec::new();
    for b in 0..BINS {
        let slice = &paired[b * per..if b + 1 == BINS { paired.len() } else { (b + 1) * per }];
        let mean_amp = slice.iter().map(|(a, _)| a).sum::<f64>() / slice.len() as f64;
        let ts: Vec<f64> = slice.iter().map(|(_, t)| *t).collect();
        let (ph, _) = forcing.phases(&ts);
        let obs = stats::schuster(&ph, 1)[0].d_squared / ph.len() as f64;

        let mut null: Vec<f64> = offsets
            .iter()
            .map(|&off| {
                let sh: Vec<f64> = ts.iter().filter_map(|&t| forcing.phase_at(t + off)).collect();
                stats::schuster(&sh, 1)[0].d_squared / sh.len() as f64
            })
            .collect();
        null.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let ge = null.iter().filter(|&&x| x >= obs).count();
        let p = (ge as f64 + 1.0) / (NULL_TRIALS as f64 + 1.0);

        rows.push((mean_amp, obs));
        println!(
            "{:>4} {:>9} {:>12.3e} {:>10.1} {:>10.1} {:>8.4}{}",
            b + 1,
            slice.len(),
            mean_amp,
            obs,
            null.last().copied().unwrap_or(0.0),
            p,
            if p < 0.05 { "  *" } else { "" }
        );
    }

    // Fit log(D2/N) against log(amplitude). Linear response predicts slope 2,
    // since D grows with amplitude and D^2 with its square.
    let (a0, d0) = rows[0];
    println!("\nscaling relative to bin 1:");
    println!("{:>4} {:>12} {:>12}", "bin", "amp ratio", "D2/N ratio");
    for (a, d) in &rows {
        println!("{:>4} {:>12.2} {:>12.2}", rows.iter().position(|r| r.0 == *a).unwrap() + 1, a / a0, d / d0);
    }

    let n_f = rows.len() as f64;
    let (sx, sy): (f64, f64) = rows.iter().fold((0.0, 0.0), |(x, y), (a, d)| {
        (x + (a / a0).ln(), y + (d / d0).ln())
    });
    let (mx, my) = (sx / n_f, sy / n_f);
    let (num, den): (f64, f64) = rows.iter().fold((0.0, 0.0), |(nu, de), (a, d)| {
        let dx = (a / a0).ln() - mx;
        (nu + dx * ((d / d0).ln() - my), de + dx * dx)
    });
    let slope = num / den;
    println!("\nlog-log slope: {slope:.2}");

    // The per-bin p-values above are NOT trustworthy, and bin 6 shows why: it has
    // the highest D2/N yet a large p. Binning by amplitude selects on the forcing
    // itself, and high-amplitude cycles recur at the ~14.77 d spring-neap period,
    // so each bin inherits temporal structure derived from the very signal under
    // test. The per-bin null does not preserve that.
    //
    // The claim is the TREND, so the trend is what must be nulled: shift the
    // event times, RE-BIN at the shifted times, and refit the slope.
    println!("\nnulling the trend itself (re-binning at shifted times):");
    // Fixed amplitude thresholds, taken once from the observed distribution and
    // reused for every null trial. Sorting 1.5M events per trial would dominate
    // runtime; fixed edges make binning O(n) and keep the binning identical
    // between observed and null, which is what a like-for-like comparison needs.
    let edges: Vec<f64> = (1..BINS)
        .map(|b| paired[b * paired.len() / BINS].0)
        .collect();
    let bin_of = |a: f64| edges.partition_point(|&e| e <= a);

    let slope_of = |ts: &[f64]| -> Option<f64> {
        let mut sum = vec![0.0f64; BINS];
        let mut cnt = vec![0usize; BINS];
        let mut re = vec![0.0f64; BINS];
        let mut im = vec![0.0f64; BINS];
        for &t in ts {
            let (Some(a), Some(ph)) = (forcing.cycle_amplitude_at(t), forcing.phase_at(t)) else {
                continue;
            };
            let b = bin_of(a).min(BINS - 1);
            sum[b] += a;
            cnt[b] += 1;
            let (s, c) = ph.sin_cos();
            re[b] += c;
            im[b] += s;
        }
        if cnt.iter().any(|&c| c < 100) {
            return None;
        }
        let rows: Vec<(f64, f64)> = (0..BINS)
            .map(|b| {
                let n = cnt[b] as f64;
                (sum[b] / n, (re[b] * re[b] + im[b] * im[b]) / n)
            })
            .collect();
        let (a0, d0) = rows[0];
        let nf = rows.len() as f64;
        let (sx, sy): (f64, f64) = rows.iter().fold((0.0, 0.0), |(x, y), (a, d)| {
            (x + (a / a0).ln(), y + (d / d0).ln())
        });
        let (mx, my) = (sx / nf, sy / nf);
        let (nu, de): (f64, f64) = rows.iter().fold((0.0, 0.0), |(nu, de), (a, d)| {
            let dx = (a / a0).ln() - mx;
            (nu + dx * ((d / d0).ln() - my), de + dx * dx)
        });
        Some(nu / de)
    };

    let mut null_slopes: Vec<f64> = offsets
        .iter()
        .filter_map(|&off| {
            let sh: Vec<f64> = times.iter().map(|&t| t + off).collect();
            slope_of(&sh)
        })
        .collect();
    null_slopes.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let ge = null_slopes.iter().filter(|&&x| x >= slope).count();
    let p_slope = (ge as f64 + 1.0) / (null_slopes.len() as f64 + 1.0);
    println!(
        "  observed slope {slope:.2}   null median {:.2}   null max {:.2}   p = {p_slope:.4}",
        null_slopes[null_slopes.len() / 2],
        null_slopes.last().copied().unwrap_or(f64::NAN)
    );
    println!(
        "  {}",
        if p_slope >= 0.05 {
            "trend not distinguishable from the null -> do not claim an amplitude law"
        } else if slope < 2.6 {
            "trend real, consistent with a linear response (slope 2)"
        } else {
            "trend real and steeper than linear -> consistent with the exponential regime"
        }
    );
    Ok(())
}
