//! A2: compute the real tidal field on the Moon and sanity-check its geometry.
//!
//!     cargo run --example lunar_tidal_field
//!
//! Requires kernels — see `scripts/fetch-kernels.sh`.
//!
//! Checks that fall out of the physics, not from any fit:
//!
//! - Earth's tidal effect on the Moon should be ~140–190× the Sun's. The nominal
//!   GM/d³ ratio at mean distances is ~178, but it swings widely: lunar distance
//!   varies 5.5% (and tides go as 1/d³), and Earth's heliocentric distance varies
//!   3.3%. Near 1 January the Moon is close to apogee *and* Earth near perihelion,
//!   both pushing the ratio down — this is the perigee–apogee amplitude knob of
//!   `docs/08-hypotheses.md` §13c showing up in the raw geometry.
//! - The principal axis should point close to Earth, since Earth dominates
//! - The tensor must be trace-free
//! - The Earth-direction angle in `MOON_PA` should stay small — the Moon is
//!   tidally locked, so Earth stays near the prime meridian, wandering only by
//!   libration (a few degrees)

use ph_core::field::TidalField;
use rustspice_core::{Aberration, Et, KernelSet};

fn main() -> rustspice_core::Result<()> {
    let mut ks = KernelSet::new();
    for k in [
        "naif0012.tls",
        "de440s.bsp",
        "pck00011.tpc",
        "gm_de440.tpc",
        "moon_pa_de440_200625.bpc",
        "moon_de440_250416.tf",
    ] {
        ks.add_file(format!("kernels/{k}"))?;
    }
    let mut spice = ks.open()?;

    let field = TidalField::on_moon(&mut spice)?;
    println!("stressing bodies:");
    for (name, gm) in field.bodies() {
        println!("  {name:<8} GM = {gm:.4} km^3/s^2");
    }

    // Mid-Apollo epoch.
    let et = spice.parse_time("1973-01-01T00:00:00")?;
    let t = field.tensor_at(&mut spice, et)?;

    println!("\ntidal tensor on the Moon at 1973-01-01, MOON_PA frame (s^-2):");
    for r in 0..3 {
        println!(
            "  [{:>12.4e} {:>12.4e} {:>12.4e}]",
            t.m[r][0], t.m[r][1], t.m[r][2]
        );
    }
    println!("  trace = {:.3e}  (analytically zero)", t.trace());

    let (vals, _) = t.eigen();
    println!("  eigenvalues = [{:.4e}, {:.4e}, {:.4e}]", vals[0], vals[1], vals[2]);

    // Relative strength of Earth versus Sun, from GM/d^3.
    let mut strength = Vec::new();
    for (name, gm) in field.bodies() {
        let p = spice.position(name, et, "MOON_PA", "MOON", Aberration::None)?;
        strength.push((name.clone(), gm / p.norm().powi(3)));
    }
    let sun = strength.iter().find(|(n, _)| n == "SUN").unwrap().1;
    println!("\nrelative tidal strength (GM/d^3):");
    for (name, s) in &strength {
        println!("  {name:<8} {:>10.1}x solar", s / sun);
    }

    // Principal axis versus the actual Earth direction.
    let axis = t.principal_axis();
    let earth = spice.position("EARTH", et, "MOON_PA", "MOON", Aberration::None)?;
    let e = [
        earth.x / earth.norm(),
        earth.y / earth.norm(),
        earth.z / earth.norm(),
    ];
    let dot: f64 = (0..3).map(|i| axis[i] * e[i]).sum();
    // The offset is the Sun: it deflects the Earth-dominated principal axis by a
    // small angle. This is the principal-axis deflection of `docs/03-tidal-tensor.md`
    // §5b, measured rather than argued.
    println!(
        "\nprincipal axis vs Earth direction: {:.3} deg off  (deflection by the Sun)",
        dot.abs().clamp(-1.0, 1.0).acos().to_degrees()
    );

    // Libration: Earth's direction in the lunar body-fixed frame over one year.
    let epochs: Vec<Et> = (0..365).map(|i| et.offset(i as f64 * 86400.0)).collect();
    let ps = spice.positions("EARTH", &epochs, "MOON_PA", "MOON", Aberration::None)?;
    let angles: Vec<f64> = ps
        .iter()
        .map(|p| (p.x / p.norm()).clamp(-1.0, 1.0).acos().to_degrees())
        .collect();
    let (lo, hi) = angles.iter().fold((f64::MAX, f64::MIN), |(l, h), &a| {
        (l.min(a), h.max(a))
    });
    println!(
        "Earth angle from lunar prime meridian over 1 yr: {lo:.2} to {hi:.2} deg  (libration)"
    );

    Ok(())
}
