//! Focal-mechanism and depth strata.
//!
//! The global scan tested every tectonic setting, depth and mechanism at once.
//! An effect confined to one of them would be diluted by the others, so the
//! obvious next move is to test them separately.
//!
//! The obvious next move is also how fishing expeditions start, and how the six
//! errors in `docs/07-research-log.md` happened. Two disciplines apply. Strata
//! are defined here, in code, before any of them is tested; and the count is kept
//! small, because splitting the sample raises the detectable effect as 1/sqrt(n)
//! and a stratum too small to detect anything is worse than useless — it adds a
//! multiplicity penalty while contributing no information.
//!
//! # Classifying a mechanism
//!
//! From the rake, which is stable under the nodal-plane ambiguity for the three
//! main classes: a thrust reads near +90 on *both* planes, a normal near −90 on
//! both, and a strike-slip near 0 or 180 on both. That is why rake is used here
//! rather than the P/T/B axis plunges of the Frohlich ternary — it needs no
//! choice between two equally valid planes.

use eqf::gcmt::Cmt;

/// Focal mechanism class.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Mechanism {
    Thrust,
    Normal,
    StrikeSlip,
}

impl Mechanism {
    pub fn name(&self) -> &'static str {
        match self {
            Mechanism::Thrust => "thrust",
            Mechanism::Normal => "normal",
            Mechanism::StrikeSlip => "strike-slip",
        }
    }
}

/// Wrap an angle to (−180, 180].
fn wrap180(x: f64) -> f64 {
    let mut r = x % 360.0;
    if r > 180.0 {
        r -= 360.0;
    }
    if r <= -180.0 {
        r += 360.0;
    }
    r
}

/// Classify a CMT by rake.
///
/// Within 45 degrees of +90 is thrust, within 45 of −90 is normal, everything
/// else is strike-slip. The 45-degree half-width splits the rake circle into
/// equal quarters, so no class is favoured by the definition.
pub fn mechanism(c: &Cmt) -> Mechanism {
    let r = wrap180(c.plane1.rake_deg);
    if (r - 90.0).abs() < 45.0 {
        Mechanism::Thrust
    } else if (r + 90.0).abs() < 45.0 {
        Mechanism::Normal
    } else {
        Mechanism::StrikeSlip
    }
}

/// Depth class. The boundary at 70 km is the conventional shallow/intermediate
/// line, and it separates crustal and interface faulting from slab processes.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Depth {
    Shallow,
    Deep,
}

pub fn depth_class(c: &Cmt) -> Depth {
    if c.depth_km < 70.0 {
        Depth::Shallow
    } else {
        Depth::Deep
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ph_core::fault::FaultPlane;

    fn cmt(rake: f64, depth: f64) -> Cmt {
        Cmt {
            day: 0.0,
            lat_deg: 0.0,
            lon_deg: 0.0,
            depth_km: depth,
            mw: 6.0,
            plane1: FaultPlane::new(0.0, 45.0, rake),
            plane2: FaultPlane::new(180.0, 45.0, rake),
        }
    }

    #[test]
    fn rake_classifies_the_three_mechanisms() {
        assert_eq!(mechanism(&cmt(90.0, 10.0)), Mechanism::Thrust);
        assert_eq!(mechanism(&cmt(70.0, 10.0)), Mechanism::Thrust);
        assert_eq!(mechanism(&cmt(-90.0, 10.0)), Mechanism::Normal);
        assert_eq!(mechanism(&cmt(-110.0, 10.0)), Mechanism::Normal);
        assert_eq!(mechanism(&cmt(0.0, 10.0)), Mechanism::StrikeSlip);
        assert_eq!(mechanism(&cmt(180.0, 10.0)), Mechanism::StrikeSlip);
        assert_eq!(mechanism(&cmt(-180.0, 10.0)), Mechanism::StrikeSlip);
    }

    #[test]
    fn the_three_classes_partition_the_rake_circle_evenly() {
        // No class may be favoured by the definition itself: over a uniform sweep
        // of rake, thrust and normal must each take a quarter and strike-slip a
        // half (it spans both 0 and 180).
        let mut counts = [0usize; 3];
        let n = 3600;
        for k in 0..n {
            let r = -180.0 + 360.0 * k as f64 / n as f64;
            match mechanism(&cmt(r, 10.0)) {
                Mechanism::Thrust => counts[0] += 1,
                Mechanism::Normal => counts[1] += 1,
                Mechanism::StrikeSlip => counts[2] += 1,
            }
        }
        let f = |i: usize| counts[i] as f64 / n as f64;
        assert!((f(0) - 0.25).abs() < 0.01, "thrust {}", f(0));
        assert!((f(1) - 0.25).abs() < 0.01, "normal {}", f(1));
        assert!((f(2) - 0.50).abs() < 0.01, "strike-slip {}", f(2));
    }

    #[test]
    fn classification_survives_the_nodal_plane_ambiguity() {
        // A pure thrust reads near +90 on both planes and a pure normal near -90
        // on both, which is why rake works without choosing a plane. Checked on
        // the conjugate pairs rather than asserted.
        for (r1, r2, want) in [
            (90.0, 90.0, Mechanism::Thrust),
            (-90.0, -90.0, Mechanism::Normal),
            (0.0, 180.0, Mechanism::StrikeSlip),
            (180.0, 0.0, Mechanism::StrikeSlip),
        ] {
            let mut c = cmt(r1, 10.0);
            c.plane2 = FaultPlane::new(0.0, 45.0, r2);
            assert_eq!(mechanism(&c), want, "rakes {r1}/{r2}");
            // Swapping which plane is "first" must not change the answer.
            let swapped = Cmt { plane1: c.plane2, plane2: c.plane1, ..c };
            assert_eq!(mechanism(&swapped), want, "swapped rakes {r1}/{r2}");
        }
    }

    #[test]
    fn depth_splits_at_the_conventional_boundary() {
        assert_eq!(depth_class(&cmt(90.0, 15.0)), Depth::Shallow);
        assert_eq!(depth_class(&cmt(90.0, 69.9)), Depth::Shallow);
        assert_eq!(depth_class(&cmt(90.0, 70.0)), Depth::Deep);
        assert_eq!(depth_class(&cmt(90.0, 400.0)), Depth::Deep);
    }
}
