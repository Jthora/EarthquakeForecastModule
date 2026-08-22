//! Declustering — removing aftershocks before anything else is asked of the data.
//!
//! # Why this cannot be skipped
//!
//! A magnitude 7 rupture is followed by hundreds of smaller events in the same
//! place over the following days. They are not independent samples of "when
//! earthquakes happen": they happened then because the mainshock happened.
//!
//! For the matched design in [`crate::sampling`] this is the worst possible
//! confounder. Controls are drawn from the same cell a few days from the case, so
//! an aftershock sequence puts hundreds of cases into a window where the *other*
//! cases' controls also sit. Any feature that varies over a few days then appears
//! to separate cases from controls — not because it drives seismicity, but because
//! the catalogue is lumpy and the lumps are about as wide as the matching window.
//!
//! That mechanism is behind a large share of the irreproducible results in the
//! tidal-triggering literature, which is why the declustered catalogue is the
//! default input here rather than an option.
//!
//! # Gardner–Knopoff
//!
//! The 1974 space-time windows, still the standard reference method:
//!
//! ```text
//! L(M) = 10^(0.1238 M + 0.983)                     km
//! T(M) = 10^(0.032 M + 2.7389)     M ≥ 6.5         days
//!      = 10^(0.5409 M − 0.547)     M < 6.5
//! ```
//!
//! Events are taken largest-magnitude first; each one not already claimed becomes a
//! mainshock and claims every smaller event inside its window. Processing by
//! magnitude rather than by time makes the result independent of catalogue order
//! and removes foreshocks as well as aftershocks.
//!
//! The windows are deliberately generous, so this **removes some independent events
//! along with the dependent ones**. That is the right direction to err: a surviving
//! aftershock manufactures signal, while a discarded independent event only costs
//! statistical power.

use eqf::comcat::Quake;

/// Gardner–Knopoff distance window, km.
pub fn window_km(magnitude: f64) -> f64 {
    10f64.powf(0.1238 * magnitude + 0.983)
}

/// Gardner–Knopoff time window, days.
pub fn window_days(magnitude: f64) -> f64 {
    if magnitude >= 6.5 {
        10f64.powf(0.032 * magnitude + 2.7389)
    } else {
        10f64.powf(0.5409 * magnitude - 0.547)
    }
}

/// Great-circle distance between two points, km.
pub fn distance_km(lat1: f64, lon1: f64, lat2: f64, lon2: f64) -> f64 {
    let (p1, p2) = (lat1.to_radians(), lat2.to_radians());
    let dp = p2 - p1;
    let dl = (lon2 - lon1).to_radians();
    // Haversine rather than the spherical law of cosines: the latter loses all
    // precision at small separations, which is exactly the regime that decides
    // whether two events belong to the same cluster.
    let a = (dp / 2.0).sin().powi(2) + p1.cos() * p2.cos() * (dl / 2.0).sin().powi(2);
    2.0 * crate::cells::EARTH_RADIUS_KM * a.sqrt().clamp(0.0, 1.0).asin()
}

/// Split a catalogue into independent events and dependent ones.
///
/// Returns `(mainshocks, dependent)`, both in the input's original order.
pub fn gardner_knopoff(quakes: &[Quake]) -> (Vec<Quake>, Vec<Quake>) {
    let n = quakes.len();
    let mut order: Vec<usize> = (0..n).collect();
    // Largest first; ties broken by time, then index, so the result never depends
    // on sort stability or on the order the catalogue happened to arrive in.
    order.sort_by(|&a, &b| {
        quakes[b]
            .magnitude
            .partial_cmp(&quakes[a].magnitude)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then(
                quakes[a]
                    .day
                    .partial_cmp(&quakes[b].day)
                    .unwrap_or(std::cmp::Ordering::Equal),
            )
            .then(a.cmp(&b))
    });

    // Time-sorted index, so each mainshock's window is found by binary search
    // instead of scanning the catalogue — the difference between minutes and days
    // on half a million events.
    let mut by_time: Vec<usize> = (0..n).collect();
    by_time.sort_by(|&a, &b| {
        quakes[a]
            .day
            .partial_cmp(&quakes[b].day)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    let times: Vec<f64> = by_time.iter().map(|&i| quakes[i].day).collect();

    let mut claimed = vec![false; n];
    let mut is_main = vec![false; n];

    for &i in &order {
        if claimed[i] {
            continue;
        }
        claimed[i] = true;
        is_main[i] = true;
        let q = &quakes[i];
        let dt = window_days(q.magnitude);
        let dr = window_km(q.magnitude);

        let lo = times.partition_point(|&t| t < q.day - dt);
        let hi = times.partition_point(|&t| t <= q.day + dt);
        for &j in &by_time[lo..hi] {
            if claimed[j] {
                continue;
            }
            if distance_km(q.lat_deg, q.lon_deg, quakes[j].lat_deg, quakes[j].lon_deg) <= dr {
                claimed[j] = true;
            }
        }
    }

    let mut main = Vec::new();
    let mut dep = Vec::new();
    for (i, q) in quakes.iter().enumerate() {
        if is_main[i] {
            main.push(*q);
        } else {
            dep.push(*q);
        }
    }
    (main, dep)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn q(day: f64, lat: f64, lon: f64, mag: f64) -> Quake {
        Quake { day, lat_deg: lat, lon_deg: lon, depth_km: 10.0, magnitude: mag }
    }

    #[test]
    fn windows_grow_with_magnitude_and_match_published_values() {
        // Spot values from Gardner & Knopoff 1974, table 1.
        assert!((window_km(6.0) - 54.0).abs() < 2.0, "{}", window_km(6.0));
        assert!((window_km(7.0) - 72.0).abs() < 3.0, "{}", window_km(7.0));
        assert!((window_days(6.0) - 510.0).abs() < 30.0, "{}", window_days(6.0));
        assert!((window_days(7.0) - 915.0).abs() < 40.0, "{}", window_days(7.0));
        for m in [4.0, 5.0, 6.0, 7.0, 8.0] {
            assert!(window_km(m) > window_km(m - 0.5));
            assert!(window_days(m) > window_days(m - 0.5));
        }
    }

    #[test]
    fn distance_is_accurate_at_both_scales() {
        // A degree of latitude at the equator.
        assert!((distance_km(0.0, 0.0, 1.0, 0.0) - 111.19).abs() < 0.05);
        // Antipodes.
        let half = std::f64::consts::PI * crate::cells::EARTH_RADIUS_KM;
        assert!((distance_km(0.0, 0.0, 0.0, 180.0) - half).abs() < 1e-6);
        // Small separations, where the law of cosines would lose all its digits.
        let d = distance_km(35.0, 139.0, 35.0, 139.001);
        assert!((d - 0.0911).abs() < 0.001, "{d}");
        assert_eq!(distance_km(12.0, 34.0, 12.0, 34.0), 0.0);
    }

    #[test]
    fn an_aftershock_sequence_collapses_to_its_mainshock() {
        let mut cat = vec![q(100.0, 35.0, 139.0, 7.2)];
        for k in 1..=200 {
            cat.push(q(100.0 + k as f64 * 0.05, 35.0 + k as f64 * 0.001, 139.0, 4.5));
        }
        let (main, dep) = gardner_knopoff(&cat);
        assert_eq!(main.len(), 1, "kept {} mainshocks", main.len());
        assert_eq!(main[0].magnitude, 7.2);
        assert_eq!(dep.len(), 200);
    }

    #[test]
    fn foreshocks_are_removed_too() {
        // Processing by magnitude rather than by time is what makes this work: a
        // small event two days before a large one in the same place is dependent,
        // and a time-ordered algorithm would keep it and discard the mainshock.
        let cat = vec![
            q(98.0, 35.0, 139.0, 5.0),
            q(100.0, 35.0, 139.0, 7.2),
            q(102.0, 35.0, 139.0, 5.1),
        ];
        let (main, dep) = gardner_knopoff(&cat);
        assert_eq!(main.len(), 1);
        assert_eq!(main[0].magnitude, 7.2);
        assert_eq!(dep.len(), 2);
    }

    #[test]
    fn independent_events_all_survive() {
        // Far apart in space and time: nothing should be removed.
        let cat: Vec<Quake> = (0..50)
            .map(|k| {
                q(
                    k as f64 * 4000.0,
                    (k % 7) as f64 * 20.0 - 60.0,
                    (k * 37 % 360) as f64 - 180.0,
                    5.0,
                )
            })
            .collect();
        let (main, dep) = gardner_knopoff(&cat);
        assert_eq!(dep.len(), 0, "removed {} independent events", dep.len());
        assert_eq!(main.len(), cat.len());
    }

    #[test]
    fn the_result_does_not_depend_on_input_order() {
        let mut cat = vec![q(100.0, 35.0, 139.0, 7.0)];
        for k in 1..40 {
            cat.push(q(100.0 + k as f64 * 0.3, 35.05, 139.02, 4.0 + (k % 5) as f64 * 0.1));
        }
        for k in 0..20 {
            cat.push(q(5000.0 + k as f64 * 300.0, -20.0, 100.0 + k as f64, 5.5));
        }
        let (a, _) = gardner_knopoff(&cat);
        let mut shuffled = cat.clone();
        shuffled.reverse();
        let (b, _) = gardner_knopoff(&shuffled);
        let key = |v: &[Quake]| {
            let mut k: Vec<String> = v
                .iter()
                .map(|x| format!("{:.6}/{:.4}", x.day, x.magnitude))
                .collect();
            k.sort();
            k
        };
        assert_eq!(key(&a), key(&b), "declustering is order-dependent");
    }

    #[test]
    fn every_event_is_either_a_mainshock_or_dependent() {
        let cat: Vec<Quake> = (0..500)
            .map(|k| {
                let f = k as f64;
                q(
                    f * 3.1 % 9000.0,
                    (f * 7.3) % 120.0 - 60.0,
                    (f * 11.7) % 360.0 - 180.0,
                    4.0 + (f * 0.37) % 3.5,
                )
            })
            .collect();
        let (main, dep) = gardner_knopoff(&cat);
        assert_eq!(main.len() + dep.len(), cat.len(), "events were lost or duplicated");
    }
}
