//! Parkfield low-frequency earthquake catalogue (Shelly, USGS).
//!
//! Shelly, D.R. (2017), *A 15 year catalog of more than 1 million low-frequency
//! earthquakes: tracking tremor and slip along the deep San Andreas Fault*,
//! *JGR: Solid Earth* 122, 3739–3753. Updated through April 2024.
//!
//! USGS-authored and therefore public domain:
//! <https://www.sciencebase.gov/catalog/item/67069991d34ef5df0d802308>
//!
//! # Why this catalogue
//!
//! **1,528,117 events in 88 families over 23 years**, every family holding at
//! least 5,333 events. Compare the Apollo deep moonquake nests, where the largest
//! held 85. Beeler & Lockner's requirement is ~10⁴ events for a robust tidal
//! correlation, so most families clear it *individually*.
//!
//! LFE **families** are the direct analogue of moonquake **nests**: repeating
//! events from one small patch, sharing a source mechanism.
//!
//! Parkfield LFEs also sit in the non-linear response regime. Thomas et al. (2012)
//! infer `Aσ₀ = 6×10⁻⁴ MPa` here, so `S_T/Aσ₀ ≈ 0.2–2` and the linearisation that
//! holds for ordinary crust does **not** apply — carry the full
//! `R = R₀exp(S_T/Aσ₀)` with `M > 1`.

use ph_core::catalog::{Catalog, Event};

/// Days from 1970-01-01 to a civil date (proleptic Gregorian).
///
/// Howard Hinnant's `days_from_civil`. Valid for any year, unlike the Apollo
/// module's hard-coded table — this catalogue spans 2001–2024 and will be reused
/// for terrestrial catalogues covering a century or more.
fn days_from_civil(y: i64, m: i64, d: i64) -> i64 {
    let y = if m <= 2 { y - 1 } else { y };
    let era = if y >= 0 { y } else { y - 399 }.div_euclid(400);
    let yoe = y - era * 400;
    let mp = if m > 2 { m - 3 } else { m + 9 };
    let doy = (153 * mp + 2) / 5 + d - 1;
    let doe = yoe * 365 + yoe / 4 - yoe / 100 + doy;
    era * 146_097 + doe - 719_468
}

/// Days from 1970-01-01 to 2000-01-01, the internal epoch.
const EPOCH_2000: i64 = 10_957;

/// One low-frequency earthquake.
#[derive(Debug, Clone, PartialEq)]
pub struct Lfe {
    /// Days since 2000-01-01T00:00 UTC.
    pub day: f64,
    /// Family identifier — the repeating-source patch.
    pub family: String,
    pub lat_deg: f64,
    pub lon_deg: f64,
    pub depth_km: f64,
    /// Summed cross-correlation across channels; a detection-strength proxy.
    pub ccsum: f64,
}

/// Parse the catalogue CSV.
///
/// Columns used: `year`, `month`, `day`, `s_of_day`, `ID`, `latitude`,
/// `longitude`, `depth`, `ccsum`. Unparseable rows are skipped.
pub fn parse_catalog(csv: &str) -> Vec<Lfe> {
    let mut lines = csv.lines();
    let header: Vec<&str> = match lines.next() {
        Some(h) => h.split(',').map(str::trim).collect(),
        None => return Vec::new(),
    };
    let col = |n: &str| header.iter().position(|h| *h == n);
    let (Some(cy), Some(cm), Some(cd), Some(cs), Some(cid), Some(clat), Some(clon), Some(cdep)) = (
        col("year"),
        col("month"),
        col("day"),
        col("s_of_day"),
        col("ID"),
        col("latitude"),
        col("longitude"),
        col("depth"),
    ) else {
        return Vec::new();
    };
    let cc = col("ccsum");
    let need = [cy, cm, cd, cs, cid, clat, clon, cdep]
        .into_iter()
        .max()
        .unwrap();

    let mut out = Vec::new();
    for line in lines {
        let f: Vec<&str> = line.split(',').map(str::trim).collect();
        if f.len() <= need {
            continue;
        }
        let (Ok(y), Ok(m), Ok(d), Ok(s)) = (
            f[cy].parse::<i64>(),
            f[cm].parse::<i64>(),
            f[cd].parse::<i64>(),
            f[cs].parse::<f64>(),
        ) else {
            continue;
        };
        let (Ok(lat), Ok(lon), Ok(dep)) = (
            f[clat].parse::<f64>(),
            f[clon].parse::<f64>(),
            f[cdep].parse::<f64>(),
        ) else {
            continue;
        };
        out.push(Lfe {
            day: (days_from_civil(y, m, d) - EPOCH_2000) as f64 + s / 86_400.0,
            family: f[cid].to_string(),
            lat_deg: lat,
            lon_deg: lon,
            depth_km: dep,
            ccsum: cc.and_then(|i| f.get(i)).and_then(|v| v.parse().ok()).unwrap_or(0.0),
        });
    }
    out
}

/// Events from one family, as a [`Catalog`] with times in days since 2000-01-01.
pub fn family(events: &[Lfe], id: &str) -> Catalog {
    let mut c = Catalog::new(format!("Parkfield LFE family {id}"));
    c.events = events
        .iter()
        .filter(|e| e.family == id)
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

/// Family identifiers with their event counts, largest first.
pub fn families(events: &[Lfe]) -> Vec<(String, usize)> {
    let mut map: std::collections::HashMap<&str, usize> = std::collections::HashMap::new();
    for e in events {
        *map.entry(e.family.as_str()).or_insert(0) += 1;
    }
    let mut v: Vec<(String, usize)> = map.into_iter().map(|(k, n)| (k.to_string(), n)).collect();
    v.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));
    v
}

/// Mean location of a family, for resolving stress onto its source patch.
pub fn family_location(events: &[Lfe], id: &str) -> Option<(f64, f64, f64)> {
    let mut n = 0.0;
    let (mut lat, mut lon, mut dep) = (0.0, 0.0, 0.0);
    for e in events.iter().filter(|e| e.family == id) {
        lat += e.lat_deg;
        lon += e.lon_deg;
        dep += e.depth_km;
        n += 1.0;
    }
    (n > 0.0).then(|| (lat / n, lon / n, dep / n))
}

#[cfg(test)]
mod tests {
    use super::*;

    const SAMPLE: &str = "\
year,month,day,s_of_day,hr,min,sec,ccsum,meancc,med_cc,seqday,ID,latitude,longitude,depth,n_chan
2001,4,6,43019.35,11,56,59.35,4.57,0.381,0.456,37,58459s,35.620,-120.185,25.75,12
2001,4,9,81061.60,22,31,1.60,4.96,0.413,0.475,40,33861s,35.530,-120.075,28.50,12
2000,1,1,0.00,0,0,0,1.0,0.1,0.1,1,33861s,35.500,-120.000,20.00,12
bad,1,1,0.00,0,0,0,1.0,0.1,0.1,1,33861s,35.5,-120.0,20.0,12
";

    #[test]
    fn civil_date_conversion_matches_known_anchors() {
        assert_eq!(days_from_civil(1970, 1, 1), 0);
        assert_eq!(days_from_civil(2000, 1, 1), EPOCH_2000);
        // 2000 was a leap year, 1900 was not.
        assert_eq!(days_from_civil(2000, 3, 1) - days_from_civil(2000, 2, 1), 29);
        assert_eq!(days_from_civil(1900, 3, 1) - days_from_civil(1900, 2, 1), 28);
        assert_eq!(days_from_civil(2024, 1, 1) - days_from_civil(2023, 1, 1), 365);
    }

    #[test]
    fn parses_rows_and_skips_bad_ones() {
        let e = parse_catalog(SAMPLE);
        assert_eq!(e.len(), 3, "{e:?}");
        // Epoch row is exactly day 0.
        let zero = e.iter().find(|x| x.day.abs() < 1e-12).unwrap();
        assert_eq!(zero.family, "33861s");
    }

    #[test]
    fn day_number_includes_seconds_of_day() {
        let e = parse_catalog(SAMPLE);
        let first = &e[0];
        let base = (days_from_civil(2001, 4, 6) - EPOCH_2000) as f64;
        assert!((first.day - (base + 43019.35 / 86400.0)).abs() < 1e-9);
    }

    #[test]
    fn groups_by_family() {
        let e = parse_catalog(SAMPLE);
        let f = families(&e);
        assert_eq!(f[0], ("33861s".to_string(), 2));
        assert_eq!(family(&e, "58459s").len(), 1);
    }

    #[test]
    fn family_location_averages_members() {
        let e = parse_catalog(SAMPLE);
        let (lat, lon, dep) = family_location(&e, "33861s").unwrap();
        assert!((lat - 35.515).abs() < 1e-9, "{lat}");
        assert!((lon + 120.0375).abs() < 1e-9, "{lon}");
        assert!((dep - 24.25).abs() < 1e-9, "{dep}");
        assert!(family_location(&e, "nope").is_none());
    }
}
