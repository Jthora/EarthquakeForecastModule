//! Cascadia tectonic tremor catalogue (PNSN / Wech).
//!
//! <https://pnsn.org/tremor> — public, no credentials. Fetch with
//! `scripts/fetch-cascadia.sh`.
//!
//! # Why this catalogue
//!
//! Every Parkfield result rests on **one location with co-located families**, so
//! C2's phase coherence was explicitly a coherence check rather than independent
//! confirmation. Cascadia is the cheapest genuine independence available:
//!
//! | | Parkfield | Cascadia |
//! |---|---|---|
//! | Setting | strike-slip transform | subduction megathrust |
//! | Geometry | vertical, right-lateral | shallow-dipping thrust |
//! | Location | 35.6 N, 120.2 W | 40–50 N, 122–125 W |
//! | Events | 1,528,117 LFEs | 678,084 tremor |
//! | Span | 2001–2024 | 2009–2024 |
//! | Detection | template matching | envelope cross-correlation |
//!
//! The detection method differs too, which matters: Parkfield's diurnal artifact
//! (D1) comes from template matching against a time-varying noise floor. A
//! different pipeline should carry a *different* artifact, so agreement between
//! the two sites is hard to explain instrumentally.
//!
//! # Time format
//!
//! The API returns RFC-2822-style stamps, `"Thu, 06 Aug 2009 00:00:00 GMT"`.
//! Parsed here to **days since 2000-01-01T00:00 UTC**, matching the other
//! catalogue modules.

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

fn month_number(name: &str) -> Option<i64> {
    const M: [&str; 12] = [
        "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec",
    ];
    M.iter().position(|m| *m == name).map(|i| i as i64 + 1)
}

/// One tremor detection.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Tremor {
    /// Days since 2000-01-01T00:00 UTC.
    pub day: f64,
    pub lat_deg: f64,
    pub lon_deg: f64,
    pub depth_km: f64,
}

/// Parse `"Thu, 06 Aug 2009 00:00:00 GMT"` to days since 2000-01-01.
pub fn parse_time(s: &str) -> Option<f64> {
    let s = s.trim().trim_matches('"');
    // Drop the weekday prefix if present.
    let rest = s.split_once(',').map(|(_, r)| r.trim()).unwrap_or(s);
    let mut it = rest.split_whitespace();
    let d: i64 = it.next()?.parse().ok()?;
    let mon = month_number(it.next()?)?;
    let y: i64 = it.next()?.parse().ok()?;
    let hms = it.next()?;
    let mut t = hms.split(':');
    let hh: f64 = t.next()?.parse().ok()?;
    let mm: f64 = t.next()?.parse().ok()?;
    let ss: f64 = t.next().unwrap_or("0").parse().ok()?;
    if !(0.0..24.0).contains(&hh) || !(0.0..60.0).contains(&mm) || !(0.0..61.0).contains(&ss) {
        return None;
    }
    Some((days_from_civil(y, mon, d) - EPOCH_2000) as f64 + (hh * 3600.0 + mm * 60.0 + ss) / 86_400.0)
}

/// Parse the CSV produced by `scripts/fetch-cascadia.sh`.
pub fn parse_catalog(csv: &str) -> Vec<Tremor> {
    let mut lines = csv.lines();
    let header: Vec<&str> = match lines.next() {
        Some(h) => h.split(',').map(str::trim).collect(),
        None => return Vec::new(),
    };
    let col = |n: &str| header.iter().position(|h| *h == n);
    let (Some(ct), Some(cla), Some(clo), Some(cd)) =
        (col("time_iso"), col("lat"), col("lon"), col("depth_km"))
    else {
        return Vec::new();
    };

    let mut out = Vec::new();
    for line in lines {
        // The time field is quoted and contains a comma, so split it off first.
        let (time_field, rest) = if let Some(stripped) = line.strip_prefix('"') {
            match stripped.split_once('"') {
                Some((t, r)) => (t, r.trim_start_matches(',')),
                None => continue,
            }
        } else {
            match line.split_once(',') {
                Some((t, r)) => (t, r),
                None => continue,
            }
        };
        let _ = ct;
        let f: Vec<&str> = rest.split(',').map(str::trim).collect();
        // Remaining columns shift left by one now that time is consumed.
        let (Some(day), Some(lat), Some(lon)) = (
            parse_time(time_field),
            f.get(cla - 1).and_then(|v| v.parse::<f64>().ok()),
            f.get(clo - 1).and_then(|v| v.parse::<f64>().ok()),
        ) else {
            continue;
        };
        out.push(Tremor {
            day,
            lat_deg: lat,
            lon_deg: lon,
            depth_km: f.get(cd - 1).and_then(|v| v.parse().ok()).unwrap_or(f64::NAN),
        });
    }
    out
}

/// Events within a latitude band, as a [`Catalog`] with times in days.
///
/// Cascadia tremor spans roughly 40–50 N. Banding by latitude gives
/// quasi-independent sub-populations along strike, the closest analogue here to
/// Parkfield's families or the Moon's nests.
pub fn latitude_band(events: &[Tremor], lat_min: f64, lat_max: f64) -> Catalog {
    let mut c = Catalog::new(format!("Cascadia tremor {lat_min:.0}-{lat_max:.0} N"));
    c.events = events
        .iter()
        .filter(|e| e.lat_deg >= lat_min && e.lat_deg < lat_max)
        .map(|e| Event {
            et: e.day,
            lat_deg: Some(e.lat_deg),
            lon_deg: Some(e.lon_deg),
            depth_km: Some(e.depth_km),
            magnitude: None,
        })
        .collect();
    c
}

/// Mean location of a latitude band.
pub fn band_location(events: &[Tremor], lat_min: f64, lat_max: f64) -> Option<(f64, f64)> {
    let (mut la, mut lo, mut n) = (0.0, 0.0, 0.0);
    for e in events.iter().filter(|e| e.lat_deg >= lat_min && e.lat_deg < lat_max) {
        la += e.lat_deg;
        lo += e.lon_deg;
        n += 1.0;
    }
    (n > 0.0).then(|| (la / n, lo / n))
}

#[cfg(test)]
mod tests {
    use super::*;

    const SAMPLE: &str = "time_iso,lat,lon,depth_km,magnitude,duration_s,energy,num_stas\n\
\"Thu, 06 Aug 2009 00:00:00 GMT\",46.63,-122.41,60.0,,,0.0,0.0\n\
\"Mon, 01 Jan 2018 00:12:30 GMT\",45.49,-122.75,20.0,1.0,150.0,41412.0,7.0\n\
\"Bad, 99 Xxx 2018 00:00:00 GMT\",45.0,-122.0,30.0,,,,\n";

    #[test]
    fn parses_rfc_style_times() {
        // 2018-01-01 is 6575 days after 2000-01-01.
        let d = parse_time("\"Mon, 01 Jan 2018 00:00:00 GMT\"").unwrap();
        assert_eq!(days_from_civil(2018, 1, 1) - EPOCH_2000, 6575);
        assert!((d - 6575.0).abs() < 1e-9, "{d}");
        // Seconds of day carry through.
        let d2 = parse_time("Mon, 01 Jan 2018 12:30:15 GMT").unwrap();
        assert!((d2 - (6575.0 + (12.0 * 3600.0 + 30.0 * 60.0 + 15.0) / 86400.0)).abs() < 1e-9);
    }

    #[test]
    fn rejects_malformed_times() {
        assert!(parse_time("nonsense").is_none());
        assert!(parse_time("Bad, 99 Xxx 2018 00:00:00 GMT").is_none());
        assert!(parse_time("Mon, 01 Jan 2018 99:00:00 GMT").is_none());
    }

    #[test]
    fn parses_catalog_and_skips_bad_rows() {
        let t = parse_catalog(SAMPLE);
        assert_eq!(t.len(), 2, "{t:?}");
        assert!((t[0].lat_deg - 46.63).abs() < 1e-9);
        assert!((t[0].lon_deg + 122.41).abs() < 1e-9);
        assert!((t[0].depth_km - 60.0).abs() < 1e-9);
    }

    #[test]
    fn bands_by_latitude() {
        let t = parse_catalog(SAMPLE);
        assert_eq!(latitude_band(&t, 46.0, 47.0).len(), 1);
        assert_eq!(latitude_band(&t, 45.0, 46.0).len(), 1);
        assert_eq!(latitude_band(&t, 48.0, 49.0).len(), 0);
        let (la, _) = band_location(&t, 46.0, 47.0).unwrap();
        assert!((la - 46.63).abs() < 1e-9);
        assert!(band_location(&t, 10.0, 11.0).is_none());
    }
}
