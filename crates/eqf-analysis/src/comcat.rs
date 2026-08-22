//! USGS ComCat global earthquake catalogue.
//!
//! <https://earthquake.usgs.gov/fdsnws/event/1/> — public, no credentials. Fetch
//! with `scripts/fetch-comcat.sh`.
//!
//! # This is the third rung of the validation ladder
//!
//! Phase 1 recovered known deep moonquake periodicities. Phase 2 measured tidal
//! response at two independent tremor sites. Both are **controls** for what this
//! catalogue is for: tremor has a short `T_a`, so its short-period response says
//! nothing about ordinary crust.
//!
//! The band prediction is that ordinary crust responds at **~1 yr to ~200 yr** and
//! is damped below. Measuring `R(ω)` here and comparing against tremor is the test.
//!
//! # ⚠ Magnitude of completeness is not a detail
//!
//! Detection capability improved enormously over any long catalogue, and that
//! trend **projects onto long-period features and manufactures signal** — precisely
//! where the band prediction looks. Choosing the threshold wrongly would fabricate
//! the result we are testing for.
//!
//! Decadal counts settle it empirically:
//!
//! | Threshold | 1970s | 1980s | 1990s | 2000s | 2010s | verdict |
//! |---|---|---|---|---|---|---|
//! | M5.0+ | 13581 | 16025 | 14660 | 17234 | 18469 | **rises 36% — incomplete** |
//! | **M5.5+** | 4377 | 4384 | 4865 | 5160 | 4898 | **18% spread, no trend — stable** |
//! | M6.0+ | 1127 | 1287 | 1535 | 1585 | 1494 | rises 33% |
//!
//! **Use M ≥ 5.5.** It gives 25,962 events over 1970–2025, comfortably above
//! Beeler & Lockner's ~13,000 requirement.
//!
//! # Longitude matters here
//!
//! Unlike the single-site catalogues, these events span the globe. Tidal phase is
//! *local*, so [`crate::doodson::Constituent::phase_at_longitude`] must be used —
//! for semidiurnal constituents a 180° longitude error is a full cycle of phase.

use ph_core::catalog::{Catalog, Event};

/// Days from 1970-01-01 to a civil date (proleptic Gregorian), Hinnant.
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

/// One catalogued earthquake.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Quake {
    /// Days since 2000-01-01T00:00 UTC.
    pub day: f64,
    pub lat_deg: f64,
    /// East longitude, degrees, as ComCat reports it (−180 to 180).
    pub lon_deg: f64,
    pub depth_km: f64,
    pub magnitude: f64,
}

/// Parse an ISO-8601 stamp like `1970-01-01T17:11:00.630Z` to days since 2000.
pub fn parse_time(s: &str) -> Option<f64> {
    let s = s.trim();
    let (date, rest) = s.split_once('T')?;
    let mut d = date.split('-');
    let y: i64 = d.next()?.parse().ok()?;
    let mo: i64 = d.next()?.parse().ok()?;
    let da: i64 = d.next()?.parse().ok()?;
    if !(1..=12).contains(&mo) || !(1..=31).contains(&da) {
        return None;
    }
    let time = rest.trim_end_matches('Z');
    let mut t = time.split(':');
    let hh: f64 = t.next()?.parse().ok()?;
    let mm: f64 = t.next()?.parse().ok()?;
    let ss: f64 = t.next().unwrap_or("0").parse().ok()?;
    // Leap seconds appear as :60 in some catalogues; allow and let them round in.
    if !(0.0..24.0).contains(&hh) || !(0.0..60.0).contains(&mm) || !(0.0..61.0).contains(&ss) {
        return None;
    }
    Some(
        (days_from_civil(y, mo, da) - EPOCH_2000) as f64
            + (hh * 3600.0 + mm * 60.0 + ss) / 86_400.0,
    )
}

/// Parse the CSV produced by `scripts/fetch-comcat.sh`.
pub fn parse_catalog(csv: &str) -> Vec<Quake> {
    let mut lines = csv.lines();
    let header: Vec<&str> = match lines.next() {
        Some(h) => h.split(',').map(str::trim).collect(),
        None => return Vec::new(),
    };
    let col = |n: &str| header.iter().position(|h| *h == n);
    let (Some(ct), Some(cla), Some(clo), Some(cd), Some(cm)) = (
        col("time"),
        col("latitude"),
        col("longitude"),
        col("depth"),
        col("mag"),
    ) else {
        return Vec::new();
    };
    let need = [ct, cla, clo, cd, cm].into_iter().max().unwrap();

    let mut out = Vec::new();
    for line in lines {
        let f: Vec<&str> = line.split(',').map(str::trim).collect();
        if f.len() <= need {
            continue;
        }
        let (Some(day), Some(lat), Some(lon)) = (
            parse_time(f[ct]),
            f[cla].parse::<f64>().ok(),
            f[clo].parse::<f64>().ok(),
        ) else {
            continue;
        };
        out.push(Quake {
            day,
            lat_deg: lat,
            lon_deg: lon,
            depth_km: f[cd].parse().unwrap_or(f64::NAN),
            magnitude: f[cm].parse().unwrap_or(f64::NAN),
        });
    }
    out
}

/// Filter by magnitude and depth, then return times sorted ascending.
pub fn select(
    quakes: &[Quake],
    min_mag: f64,
    max_depth_km: Option<f64>,
) -> Vec<Quake> {
    let mut v: Vec<Quake> = quakes
        .iter()
        .copied()
        .filter(|q| q.magnitude >= min_mag)
        .filter(|q| max_depth_km.is_none_or(|d| q.depth_km <= d))
        .collect();
    v.sort_by(|a, b| a.day.partial_cmp(&b.day).unwrap());
    v
}

/// Convert to a [`Catalog`], with times in days since 2000-01-01.
pub fn to_catalog(quakes: &[Quake], name: impl Into<String>) -> Catalog {
    let mut c = Catalog::new(name);
    c.events = quakes
        .iter()
        .map(|q| Event {
            et: q.day,
            lat_deg: Some(q.lat_deg),
            lon_deg: Some(q.lon_deg),
            depth_km: Some(q.depth_km),
            magnitude: Some(q.magnitude),
        })
        .collect();
    c
}

/// Event counts per decade at a magnitude threshold — the completeness check.
///
/// A threshold is usable when counts show **no monotonic trend**. A rising series
/// means the early catalogue is incomplete, and that incompleteness will appear as
/// long-period signal.
pub fn decadal_counts(quakes: &[Quake], min_mag: f64, start_year: i64, decades: usize) -> Vec<usize> {
    let d0 = (days_from_civil(start_year, 1, 1) - EPOCH_2000) as f64;
    (0..decades)
        .map(|i| {
            let lo = d0 + i as f64 * 3652.5;
            let hi = lo + 3652.5;
            quakes
                .iter()
                .filter(|q| q.magnitude >= min_mag && q.day >= lo && q.day < hi)
                .count()
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    const SAMPLE: &str = "time,latitude,longitude,depth,mag\n\
1970-01-01T17:11:00.630Z,-29.4,-177.169,35,5.61\n\
2000-01-01T00:00:00.000Z,0.0,0.0,10,6.0\n\
2010-06-15T12:30:00.000Z,35.5,-120.1,8,5.5\n\
not-a-time,1,2,3,4\n\
2010-13-45T00:00:00.000Z,1,2,3,4\n";

    #[test]
    fn parses_iso_times_against_known_anchors() {
        assert!((parse_time("2000-01-01T00:00:00.000Z").unwrap()).abs() < 1e-9);
        // 2010-06-15 is 3818 days after 2000-01-01.
        let d = parse_time("2010-06-15T12:30:00.000Z").unwrap();
        assert!((d - (3818.0 + 12.5 / 24.0)).abs() < 1e-9, "{d}");
        // 1970 is before the epoch, so negative.
        assert!(parse_time("1970-01-01T00:00:00.000Z").unwrap() < -10_000.0);
    }

    #[test]
    fn rejects_malformed_times() {
        assert!(parse_time("not-a-time").is_none());
        assert!(parse_time("2010-13-45T00:00:00.000Z").is_none());
        assert!(parse_time("2010-01-01T99:00:00Z").is_none());
    }

    #[test]
    fn parses_catalog_and_skips_bad_rows() {
        let q = parse_catalog(SAMPLE);
        assert_eq!(q.len(), 3, "{q:?}");
        assert!((q[0].lat_deg + 29.4).abs() < 1e-9);
        assert!((q[0].magnitude - 5.61).abs() < 1e-9);
    }

    #[test]
    fn select_filters_and_sorts() {
        let q = parse_catalog(SAMPLE);
        let s = select(&q, 5.6, None);
        assert_eq!(s.len(), 2);
        assert!(s[0].day < s[1].day, "must be time-sorted");
        // Depth cut drops the 35 km event.
        assert_eq!(select(&q, 0.0, Some(20.0)).len(), 2);
    }

    #[test]
    fn decadal_counts_bin_by_decade() {
        let q = parse_catalog(SAMPLE);
        // 2000s decade contains the 2000 event; 2010s contains the 2010 one.
        let c = decadal_counts(&q, 0.0, 2000, 2);
        assert_eq!(c, vec![1, 1]);
    }

    #[test]
    fn to_catalog_carries_location_and_magnitude() {
        let c = to_catalog(&parse_catalog(SAMPLE), "test");
        assert_eq!(c.len(), 3);
        assert!(c.events[0].lon_deg.is_some());
        assert!(c.events[0].magnitude.is_some());
    }
}
