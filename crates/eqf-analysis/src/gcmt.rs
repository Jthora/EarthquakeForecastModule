//! Global CMT focal mechanisms.
//!
//! <https://www.globalcmt.org/> — public, no credentials. Fetch with
//! `scripts/fetch-gcmt.sh`.
//!
//! # Why mechanisms matter here
//!
//! P3.4 and P3.5 tested ordinary crust using raw tidal phase and found nothing,
//! bounding the response below ~4%. But **raw tidal phase is blind to whether the
//! tide loads or unloads each fault.** A tide that promotes failure on a thrust in
//! Japan simultaneously inhibits it on a normal fault in Greece; pooled over a
//! global catalogue, those responses cancel.
//!
//! Parkfield's positive result used ΔCFS resolved on a *known* fault plane. Giving
//! the global catalogue the same treatment raises signal by **aligning the feature
//! with the physics**, not by discarding events — the opposite of the
//! stratification that P3.5 showed makes bounds worse.
//!
//! # NDK format
//!
//! Five lines per event. This parser uses:
//!
//! | Line | Fields taken |
//! |---|---|
//! | 1 | reference date, time, hypocentre |
//! | 3 | centroid time shift, centroid lat/lon/depth |
//! | 4 | moment-tensor exponent |
//! | 5 | scalar moment, **both nodal planes** (strike, dip, rake) |
//!
//! # ⚠ Nodal plane ambiguity is real and unresolvable from the mechanism alone
//!
//! A moment tensor gives **two** planes that fit equally well. Only one is the
//! fault; nothing in the CMT solution says which. This parser returns both, and
//! analyses should run each separately as a robustness check rather than silently
//! picking one. Using the wrong plane for roughly half the events dilutes a real
//! signal but does not manufacture one.

use ph_core::fault::FaultPlane;

fn days_from_civil(y: i64, m: i64, d: i64) -> i64 {
    let y = if m <= 2 { y - 1 } else { y };
    let era = if y >= 0 { y } else { y - 399 }.div_euclid(400);
    let yoe = y - era * 400;
    let mp = if m > 2 { m - 3 } else { m + 9 };
    let doy = (153 * mp + 2) / 5 + d - 1;
    let doe = yoe * 365 + yoe / 4 - yoe / 100 + doy;
    era * 146_097 + doe - 719_468
}

const EPOCH_2000: i64 = 10_957;

/// One centroid moment tensor solution.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Cmt {
    /// Centroid time, days since 2000-01-01T00:00 UTC.
    pub day: f64,
    /// Centroid latitude, degrees.
    pub lat_deg: f64,
    /// Centroid east longitude, degrees.
    pub lon_deg: f64,
    /// Centroid depth, km.
    pub depth_km: f64,
    /// Moment magnitude, from the scalar moment.
    pub mw: f64,
    /// First nodal plane.
    pub plane1: FaultPlane,
    /// Second nodal plane. Equally consistent with the mechanism.
    pub plane2: FaultPlane,
}

impl Cmt {
    /// Both nodal planes, for running an analysis twice.
    pub fn planes(&self) -> [FaultPlane; 2] {
        [self.plane1, self.plane2]
    }
}

/// Parse an NDK file.
///
/// Malformed or truncated five-line groups are skipped rather than failing the
/// parse — the concatenated monthly files occasionally carry stray blank lines.
pub fn parse_ndk(text: &str) -> Vec<Cmt> {
    let lines: Vec<&str> = text.lines().filter(|l| !l.trim().is_empty()).collect();
    let mut out = Vec::new();

    for chunk in lines.chunks(5) {
        if chunk.len() < 5 {
            break;
        }
        let Some(c) = parse_event(chunk) else { continue };
        out.push(c);
    }
    out
}

fn parse_event(l: &[&str]) -> Option<Cmt> {
    // Line 1: "PDEW 2024/01/01 07:10:09.5  37.49  137.27  10.0 0.0 7.5 NEAR ..."
    let f1: Vec<&str> = l[0].split_whitespace().collect();
    if f1.len() < 6 {
        return None;
    }
    let mut date = f1[1].split('/');
    let y: i64 = date.next()?.parse().ok()?;
    let mo: i64 = date.next()?.parse().ok()?;
    let da: i64 = date.next()?.parse().ok()?;
    if !(1..=12).contains(&mo) || !(1..=31).contains(&da) {
        return None;
    }
    let mut hms = f1[2].split(':');
    let hh: f64 = hms.next()?.parse().ok()?;
    let mm: f64 = hms.next()?.parse().ok()?;
    let ss: f64 = hms.next()?.parse().ok()?;
    let ref_day =
        (days_from_civil(y, mo, da) - EPOCH_2000) as f64 + (hh * 3600.0 + mm * 60.0 + ss) / 86_400.0;

    // Line 3: "CENTROID:  29.5 0.1  37.49 0.00  137.16 0.00  12.0  0.0 FIX ..."
    let f3: Vec<&str> = l[2].split_whitespace().collect();
    if f3.len() < 8 || !f3[0].starts_with("CENTROID") {
        return None;
    }
    let shift: f64 = f3[1].parse().ok()?;
    let lat: f64 = f3[3].parse().ok()?;
    let lon: f64 = f3[5].parse().ok()?;
    let depth: f64 = f3[7].parse().ok()?;

    // Line 4 begins with the moment-tensor exponent.
    let exponent: f64 = l[3].split_whitespace().next()?.parse().ok()?;

    // Line 5: version, three eigen triples, scalar moment, then both planes.
    let f5: Vec<&str> = l[4].split_whitespace().collect();
    if f5.len() < 17 {
        return None;
    }
    let m0_mant: f64 = f5[10].parse().ok()?;
    let p = |i: usize| -> Option<f64> { f5[i].parse().ok() };
    let plane1 = FaultPlane::new(p(11)?, p(12)?, p(13)?);
    let plane2 = FaultPlane::new(p(14)?, p(15)?, p(16)?);

    // M0 in dyne-cm; Mw = (2/3)(log10 M0 - 16.1).
    let m0 = m0_mant * 10f64.powf(exponent);
    if m0 <= 0.0 {
        return None;
    }
    let mw = (2.0 / 3.0) * (m0.log10() - 16.1);

    Some(Cmt {
        day: ref_day + shift / 86_400.0,
        lat_deg: lat,
        lon_deg: lon,
        depth_km: depth,
        mw,
        plane1,
        plane2,
    })
}

/// Filter by magnitude and depth, sorted by time.
pub fn select(cmts: &[Cmt], min_mw: f64, max_depth_km: Option<f64>) -> Vec<Cmt> {
    let mut v: Vec<Cmt> = cmts
        .iter()
        .copied()
        .filter(|c| c.mw >= min_mw)
        .filter(|c| max_depth_km.is_none_or(|d| c.depth_km <= d))
        .collect();
    v.sort_by(|a, b| a.day.partial_cmp(&b.day).unwrap());
    v
}

#[cfg(test)]
mod tests {
    use super::*;

    const SAMPLE: &str = "\
PDEW 2024/01/01 07:10:09.5  37.49  137.27  10.0 0.0 7.5 NEAR WEST COAST OF HONSH
C202401010710A   B:169  463  50 S:170  471  50 M:169  474 150 CMT: 1 TRIHD: 14.2
CENTROID:     29.5 0.1  37.49 0.00  137.16 0.00  12.0  0.0 FIX  S-20240330203700
27  2.260 0.004 -0.755 0.004 -1.500 0.004 -0.436 0.044 -0.441 0.041 -1.300 0.003
V10   2.342 82 144   0.222  2  38  -2.558  7 307   2.450  35 38   86 219 52   93
";

    #[test]
    fn parses_a_known_event() {
        let c = parse_ndk(SAMPLE);
        assert_eq!(c.len(), 1);
        let e = c[0];
        // Centroid location, not hypocentre: longitude 137.16 rather than 137.27.
        assert!((e.lat_deg - 37.49).abs() < 1e-9);
        assert!((e.lon_deg - 137.16).abs() < 1e-9, "{}", e.lon_deg);
        assert!((e.depth_km - 12.0).abs() < 1e-9);
    }

    #[test]
    fn recovers_both_nodal_planes() {
        let e = parse_ndk(SAMPLE)[0];
        assert_eq!(e.plane1, FaultPlane::new(35.0, 38.0, 86.0));
        assert_eq!(e.plane2, FaultPlane::new(219.0, 52.0, 93.0));
        // Rake near +90 on both: this is the 2024 Noto thrust event.
        assert!(e.plane1.rake_deg > 60.0 && e.plane2.rake_deg > 60.0);
        assert_eq!(e.planes().len(), 2);
    }

    #[test]
    fn magnitude_matches_the_published_value() {
        // M0 = 2.450e27 dyne-cm -> Mw = (2/3)(log10(M0) - 16.1) = 7.5
        let e = parse_ndk(SAMPLE)[0];
        assert!((e.mw - 7.5).abs() < 0.05, "Mw came out {}", e.mw);
    }

    #[test]
    fn centroid_time_includes_the_shift() {
        let e = parse_ndk(SAMPLE)[0];
        let base = (days_from_civil(2024, 1, 1) - EPOCH_2000) as f64
            + (7.0 * 3600.0 + 10.0 * 60.0 + 9.5) / 86400.0;
        assert!((e.day - (base + 29.5 / 86400.0)).abs() < 1e-9);
    }

    #[test]
    fn skips_truncated_and_malformed_groups() {
        assert!(parse_ndk("PDEW 2024/01/01 07:10:09.5  37.49  137.27  10.0").is_empty());
        let bad = SAMPLE.replace("CENTROID:", "NOTACENT:");
        assert!(parse_ndk(&bad).is_empty());
    }

    #[test]
    fn select_filters_and_sorts() {
        let c = parse_ndk(SAMPLE);
        assert_eq!(select(&c, 7.0, None).len(), 1);
        assert_eq!(select(&c, 8.0, None).len(), 0);
        assert_eq!(select(&c, 0.0, Some(5.0)).len(), 0);
    }
}
