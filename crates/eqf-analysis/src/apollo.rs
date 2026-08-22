//! Apollo Passive Seismic Experiment event catalogue.
//!
//! Parses the expanded event catalogue (`levent.1008weber.csv`) from the PDS
//! Geosciences Node, compiled by Renee Weber:
//!
//! <https://pds-geosciences.wustl.edu/lunar/urn-nasa-pds-apollo_seismic_event_catalog/data>
//!
//! # Classification: use `T2`, not `T1`
//!
//! The catalogue carries two classification columns. **`T1` is the original
//! classification and `T2` the revised one** (Nakamura 2005). They disagree for
//! 1,471 events — many originally logged as meteoroid impacts (`M`) were
//! reclassified as deep moonquakes (`A`).
//!
//! Counting `T1` alone finds 1,359 deep moonquakes; `T2` finds **7,082**, which is
//! the figure the literature reports. Always prefer `T2`.
//!
//! | Code | Event type |
//! |---|---|
//! | `A` | Deep moonquake (with nest number in `N2`) |
//! | `C` | Unclassified / ambiguous |
//! | `M` | Meteoroid impact |
//! | `Z` | Shallow moonquake |
//! | `T`, `H`, `X`, `S`, `L` | Other, including artificial impacts |
//!
//! # Time
//!
//! Times are `Y` (two-digit year, 1969–1977), `JD` (day of year, 1-based) and `S`
//! (`HHMM`, UTC). This module yields **days since 1969-01-01T00:00 UTC**, which is
//! sufficient for periodicity analysis. Conversion to ephemeris time for tidal
//! phase work is a separate step and needs a leap-second kernel.

use ph_core::catalog::{Catalog, Event};

/// Days from 1969-01-01 to 1 January of `year` (four-digit, 1969–1977).
fn days_to_year_start(year: i32) -> Option<f64> {
    const OFFSETS: [(i32, f64); 9] = [
        (1969, 0.0),
        (1970, 365.0),
        (1971, 730.0),
        (1972, 1095.0),
        (1973, 1461.0), // 1972 was a leap year
        (1974, 1826.0),
        (1975, 2191.0),
        (1976, 2556.0),
        (1977, 2922.0), // 1976 was a leap year
    ];
    OFFSETS.iter().find(|(y, _)| *y == year).map(|(_, d)| *d)
}

/// One parsed catalogue row.
#[derive(Debug, Clone, PartialEq)]
pub struct ApolloEvent {
    /// Days since 1969-01-01T00:00 UTC.
    pub day: f64,
    /// Revised classification (`T2`).
    pub class: String,
    /// Nest number (`N2`), where the event is a located deep moonquake.
    pub nest: Option<u32>,
}

/// Parse the expanded event catalogue.
///
/// Rows with an unparseable year, day-of-year or time are skipped rather than
/// failing the whole parse — the catalogue has incomplete rows by design.
pub fn parse_levent(csv: &str) -> Vec<ApolloEvent> {
    let mut lines = csv.lines();
    let header: Vec<&str> = match lines.next() {
        Some(h) => h.split(',').map(|s| s.trim()).collect(),
        None => return Vec::new(),
    };
    let col = |name: &str| header.iter().position(|h| *h == name);
    let (ci_y, ci_jd, ci_s, ci_t2, ci_n2) =
        match (col("Y"), col("JD"), col("S"), col("T2"), col("N2")) {
            (Some(a), Some(b), Some(c), Some(d), Some(e)) => (a, b, c, d, e),
            _ => return Vec::new(),
        };

    let mut out = Vec::new();
    for line in lines {
        let f: Vec<&str> = line.split(',').map(|s| s.trim()).collect();
        if f.len() <= ci_n2.max(ci_t2) {
            continue;
        }
        let (yy, jd, hhmm) = match (
            f[ci_y].parse::<i32>(),
            f[ci_jd].parse::<f64>(),
            f[ci_s].parse::<u32>(),
        ) {
            (Ok(a), Ok(b), Ok(c)) => (a, b, c),
            _ => continue,
        };
        let year = 1900 + yy;
        let Some(base) = days_to_year_start(year) else {
            continue;
        };
        let (hh, mm) = (hhmm / 100, hhmm % 100);
        if hh > 23 || mm > 59 {
            continue;
        }
        let day = base + (jd - 1.0) + (hh as f64 + mm as f64 / 60.0) / 24.0;

        out.push(ApolloEvent {
            day,
            class: f[ci_t2].to_string(),
            nest: f[ci_n2].parse::<u32>().ok(),
        });
    }
    out
}

/// Deep moonquakes (`T2 == "A"`), optionally restricted to one nest.
pub fn deep_moonquakes(events: &[ApolloEvent], nest: Option<u32>) -> Catalog {
    let name = match nest {
        Some(n) => format!("Apollo deep moonquakes, nest A{n}"),
        None => "Apollo deep moonquakes".to_string(),
    };
    let mut c = Catalog::new(name);
    c.events = events
        .iter()
        .filter(|e| e.class == "A")
        .filter(|e| nest.is_none() || e.nest == nest)
        .map(|e| Event {
            // Days, not ET. Periodicity analysis only needs relative times; see
            // the module note on time.
            et: e.day,
            lat_deg: None,
            lon_deg: None,
            depth_km: None,
            magnitude: None,
        })
        .collect();
    c
}

#[cfg(test)]
mod tests {
    use super::*;

    const SAMPLE: &str = "\
Y,JD,S,E,A1,T1,N1,T2,N2
69,208,2348,0040,3.1,M,,M,
70,001,0000,0100,1.0,M,,A,1
72,060,1200,1300,2.0,A,1,A,1
77,100,0630,0700,2.0,A,7,A,7
69,999,9999,0000,1.0,A,1,A,1
99,001,0000,0100,1.0,A,1,A,1
";

    #[test]
    fn parses_valid_rows_and_skips_bad_ones() {
        let e = parse_levent(SAMPLE);
        // Row 5 has an impossible time (9999), row 6 an out-of-range year.
        assert_eq!(e.len(), 4, "{e:?}");
    }

    #[test]
    fn day_number_accounts_for_leap_years() {
        let e = parse_levent(SAMPLE);
        // 1970-001 00:00 is exactly 365 days after 1969-01-01.
        let y70 = e.iter().find(|x| x.nest == Some(1) && x.day > 364.0).unwrap();
        assert!((y70.day - 365.0).abs() < 1e-9, "{}", y70.day);
        // 1977-100 06:30 → 2922 + 99 + 6.5/24
        let y77 = e.iter().find(|x| x.nest == Some(7)).unwrap();
        assert!((y77.day - (2922.0 + 99.0 + 6.5 / 24.0)).abs() < 1e-9);
    }

    #[test]
    fn uses_revised_classification() {
        let e = parse_levent(SAMPLE);
        // Row 2 is T1=M but T2=A — it must count as a deep moonquake.
        let dm = deep_moonquakes(&e, None);
        assert_eq!(dm.len(), 3);
    }

    #[test]
    fn filters_by_nest() {
        let e = parse_levent(SAMPLE);
        assert_eq!(deep_moonquakes(&e, Some(1)).len(), 2);
        assert_eq!(deep_moonquakes(&e, Some(7)).len(), 1);
    }
}

/// A deep moonquake nest location from `nakamura_2005_dm_locations.csv`.
///
/// Depths run 700–1200 km, but **depth does not affect the degree-2 tidal
/// tensor**: the tide-generating potential goes as `r²·P₂(cos ψ)`, so its second
/// derivatives are *constant throughout the body*. Only latitude and longitude
/// matter, and only because they set the local frame orientation. (Depth does
/// affect the elastic response — the Love-number scale factor — but that is a
/// separate correction and does not change timing.)
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct NestLocation {
    pub nest: u32,
    pub lat_deg: f64,
    pub lon_deg: f64,
    pub depth_km: f64,
    /// True where the published location is assumed rather than determined.
    pub assumed: bool,
}

/// Parse `nakamura_2005_dm_locations.csv`.
pub fn parse_dm_locations(csv: &str) -> Vec<NestLocation> {
    let mut lines = csv.lines();
    let header: Vec<&str> = match lines.next() {
        Some(h) => h.split(',').map(|s| s.trim()).collect(),
        None => return Vec::new(),
    };
    let col = |n: &str| header.iter().position(|h| *h == n);
    let (ci_a, ci_lat, ci_lon, ci_d, ci_as) = match (
        col("A"),
        col("Lat"),
        col("Long"),
        col("Depth"),
        col("Assumed"),
    ) {
        (Some(a), Some(b), Some(c), Some(d), Some(e)) => (a, b, c, d, e),
        _ => return Vec::new(),
    };

    let mut out = Vec::new();
    for line in lines {
        let f: Vec<&str> = line.split(',').map(|s| s.trim()).collect();
        if f.len() <= ci_as {
            continue;
        }
        let (Ok(nest), Ok(lat), Ok(lon), Ok(depth)) = (
            f[ci_a].parse::<u32>(),
            f[ci_lat].parse::<f64>(),
            f[ci_lon].parse::<f64>(),
            f[ci_d].parse::<f64>(),
        ) else {
            continue;
        };
        out.push(NestLocation {
            nest,
            lat_deg: lat,
            lon_deg: lon,
            depth_km: depth,
            assumed: f[ci_as].eq_ignore_ascii_case("Y"),
        });
    }
    out
}

#[cfg(test)]
mod location_tests {
    use super::*;

    const SAMPLE: &str = "\
A,Side,Lat,Lat_err,Long,Long_err,Depth,Depth_err,Assumed
1,N,-15.7,2.4,-36.6,4.6,867,29,N
5,N,1.1,94.2,-44.7,16.4,933,109,Y
bad,N,0,0,0,0,0,0,N
";

    #[test]
    fn parses_locations_and_assumed_flag() {
        let l = parse_dm_locations(SAMPLE);
        assert_eq!(l.len(), 2);
        assert_eq!(l[0].nest, 1);
        assert!((l[0].lat_deg + 15.7).abs() < 1e-9);
        assert!((l[0].lon_deg + 36.6).abs() < 1e-9);
        assert!(!l[0].assumed);
        assert!(l[1].assumed);
    }
}
