//! F2 applied: what does C4's measured band limit imply for `T_a`?
//!
//!     cargo run --release --example parkfield_ta_estimate
//!
//! # ⚠ Superseded — the premise was wrong
//!
//! This was written when C4's raw `D²/N` looked band-limited: 17,975 at M2 against
//! 338 at Mf, a 53× decline read as the falling limb above `T_a`, bounding
//! `T_a ≲ 1 d`.
//!
//! **Normalising by the forcing amplitude removed that reading.** The tidal
//! potential's own amplitudes fall with period, and once `D²/N` is divided by the
//! ΔCFS amplitude at each constituent, response per unit stress is flat to within
//! a factor of ~3 from 0.5 d to 27 d. There is no measured band limit, so there is
//! no `T_a` bound, and the "tension" this example reports **does not exist**.
//!
//! Kept because the calibration numbers it prints are correct and independently
//! useful — in particular that the M2 solid Earth tide comes out **595 Pa**
//! against Thomas et al.'s `Aσ₀ = 600 Pa`, matching to 1%, which is why Parkfield
//! sits in the non-linear regime (`S_T/Aσ₀ ≈ 1`) and why C3 measured a slope
//! steeper than linear.
//!
//! The `T_a` arithmetic below is retained as an illustration of the relation, not
//! as a result. **Do not quote the tension paragraph.**

use ph_core::love::{critical_period, stressing_rate_from, Elastic};

const YEAR_D: f64 = 365.25;

fn main() {
    let e = Elastic::EARTH;
    println!("Earth elastic calibration");
    println!("  strain factor 2h2-6l2 = {:.4}", e.strain_factor());
    println!("  stress per unit tensor = {:.3e} Pa per s^-2", e.stress_per_tensor());

    // Moon-on-Earth tensor scale, GM/d^3.
    let t_m2 = 4.9028e12 / 3.844e8_f64.powi(3);
    println!("\nM2 solid Earth tide");
    println!("  tensor    {t_m2:.3e} s^-2");
    println!("  strain    {:.3e}", e.strain(t_m2));
    println!("  stress    {:.0} Pa  ({:.2e} MPa)", e.stress(t_m2), e.stress(t_m2) / 1e6);

    // C4's bound: response still falling at 0.5 d, so the peak is at or below it.
    let t_a_days = 1.0;
    let t_a_years = t_a_days / YEAR_D;
    println!("\nC4 bound: T_a <~ {t_a_days} d");

    println!("\nimplied stressing rate, given published A*sigma0:");
    for (label, a_sigma_mpa) in [
        ("Thomas et al. 2012, Parkfield tremor", 6.0e-4),
        ("aftershock studies, low end", 1.0e-2),
        ("aftershock studies, high end", 1.0e-1),
    ] {
        let rate = stressing_rate_from(a_sigma_mpa * 1e6, t_a_years);
        println!(
            "  {label:<38} A*sigma0 {a_sigma_mpa:.0e} MPa -> tau-dot {:.2e} Pa/yr ({:.1e} MPa/yr)",
            rate,
            rate / 1e6
        );
    }

    println!("\nimplied A*sigma0, given plausible stressing rates:");
    for (label, rate_pa_yr) in [
        ("secular tectonic, shallow crust", 3.0e3),
        ("deep SAF creep loading", 1.0e4),
        ("during a slow slip episode", 1.0e5),
    ] {
        let a_sigma = rate_pa_yr * t_a_years / std::f64::consts::TAU;
        println!(
            "  {label:<38} tau-dot {rate_pa_yr:.0e} Pa/yr -> A*sigma0 {:.2e} Pa ({:.1e} MPa)",
            a_sigma,
            a_sigma / 1e6
        );
    }

    println!("\nfor contrast, ordinary crust:");
    for a in [1.0e4, 1.0e5] {
        let t = critical_period(a, 3.0e3);
        println!(
            "  A*sigma0 {:.0e} Pa, tau-dot 3e3 Pa/yr -> T_a {:.0} yr",
            a, t
        );
    }

    println!(
        "\nTENSION: T_a <~ 1 d with Thomas et al.'s A*sigma0 needs tau-dot ~1.4 MPa/yr,\n\
         some 460x secular tectonic loading. Inverting instead, a plausible tau-dot\n\
         implies A*sigma0 two orders below the published Parkfield value.\n\
         Either the band-limit inference is wrong, the two-regime model does not\n\
         describe LFEs, or one of the published numbers does not apply here.\n\
         The stress calibration is good to ~2x and cannot absorb 100x."
    );
}
