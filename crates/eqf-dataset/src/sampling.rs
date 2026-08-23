//! Case-control sampling — constructing the negative class.
//!
//! # Why not "cells with no earthquakes"
//!
//! The obvious negative class is places that never rupture. It is the wrong one.
//! A cell in the middle of a craton contributes no information about *timing*: it
//! is quiet in 1983 and quiet in 2011 for reasons that have nothing to do with any
//! feature computed here. Including it teaches a model only to recognise cratons,
//! which is the background rate restated, and its accuracy would look excellent
//! while its forecast skill stayed at zero.
//!
//! The informative negative is a place that *does* rupture, at a time when it
//! didn't. So a control is the same cell as its case, at a different instant.
//!
//! # What each scheme conditions out
//!
//! Because controls share the case's cell, everything spatial cancels: tectonic
//! setting, station density, regional magnitude bias. The offset scheme decides
//! what else cancels:
//!
//! | scheme | also conditions out | so it cannot detect |
//! |---|---|---|
//! | [`Scheme::DayOffset`] | time of day, season, network trend | anything diurnal or slower than the offset span |
//! | [`Scheme::Window`] | season, network trend | anything slower than the window |
//! | [`Scheme::Uniform`] | nothing beyond the cell | — but is exposed to every temporal bias |
//!
//! `DayOffset` is the strictest. Whole-day offsets hold local solar time fixed, so
//! the daily cycle in detection threshold — quarry blasts, traffic, cultural noise
//! — is identical for case and control and cannot be mistaken for signal. The Moon
//! meanwhile advances 12.2° per day, so lunar features still vary. The cost is
//! stated plainly in the table: this design is blind to a genuine diurnal effect,
//! and blind to slow planetary configurations that barely move in a few days.
//!
//! No single scheme is safe on its own. Running two and requiring agreement is the
//! point — a result that appears under `Uniform` but vanishes under `DayOffset` is
//! a temporal-bias artefact, and the design says so before the data is seen.

use crate::cells::Grid;
use eqf::comcat::Quake;

/// Civil date from days since 2000-01-01, by Howard Hinnant's algorithm.
///
/// Needed because the time-stratified scheme's referent window is a *calendar*
/// month, which is what makes its exchangeability argument work.
pub fn civil_from_days(z: i64) -> (i64, u32, u32) {
    let z = z + 730_425; // shift epoch to 0000-03-01
    let era = if z >= 0 { z } else { z - 146_096 } / 146_097;
    let doe = (z - era * 146_097) as i64;
    let yoe = (doe - doe / 1460 + doe / 36_524 - doe / 146_096) / 365;
    let y = yoe + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = (doy - (153 * mp + 2) / 5 + 1) as u32;
    let m = if mp < 10 { mp + 3 } else { mp - 9 } as u32;
    (if m <= 2 { y + 1 } else { y }, m, d)
}

/// How control times are drawn relative to their case.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Scheme {
    /// Offsets of ±1..=`max_days` whole days. Holds local solar time fixed.
    DayOffset { max_days: u32 },
    /// Uniform within ±`window_days` of the case.
    Window { window_days: f64 },
    /// Uniform across the whole study span.
    Uniform,
    /// Time-stratified: every time in the case's own calendar month that differs
    /// from it by a whole multiple of `spacing_days`.
    ///
    /// This is the only scheme here whose null is exactly right, and the reason is
    /// worth stating. The others place the case at the centre of its controls: with
    /// offsets of plus-or-minus k days the controls straddle the case, so for any
    /// smoothly-varying feature the case sits near its own stratum's mean while a
    /// control does not. The case is then identifiable from the feature values
    /// alone, the rows are not exchangeable, and a permutation test that assumes
    /// they are is mis-calibrated -- conservative here rather than anti-conservative,
    /// but wrong either way. Measured on the real M4 catalogue, the observed z field
    /// came out at sd 0.63 where the permutation null said 1.001 +/- 0.031.
    ///
    /// Here the referent set is fixed by the calendar, not by the case: it is every
    /// day of that month congruent to the case modulo the spacing. Under a uniform
    /// hazard the event is equally likely to have fallen on any of them, so the case
    /// really is exchangeable with its controls and the permutation null is exact.
    /// This is the standard time-stratified case-crossover referent scheme, adopted
    /// precisely because the symmetric bidirectional design has the bias above.
    ///
    /// A spacing of 7 also holds day of week fixed, so the weekly cycle in cultural
    /// noise -- and hence in detection threshold -- cannot masquerade as signal.
    /// Control count varies between 3 and 4 depending on the month.
    TimeStratified { spacing_days: u32 },
}

/// One row of the training set.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Row {
    /// Days since 2000-01-01T00:00 UTC.
    pub day: f64,
    pub cell: usize,
    /// Cell centre latitude, degrees — where site-local features are evaluated.
    pub lat_deg: f64,
    pub lon_deg: f64,
    /// True for a case (an earthquake), false for a control.
    pub case: bool,
    /// Index of the case this row belongs to. Cases and their controls share it,
    /// which is what lets a conditional model stratify on the matched set.
    pub stratum: u32,
    /// Magnitude for a case; `f64::NAN` for a control.
    pub magnitude: f64,
}

/// xoshiro256** — a small, fast, well-tested generator.
///
/// Carried explicitly rather than taken from the environment so that a run is
/// reproducible from its seed alone. A dataset that cannot be rebuilt bit-for-bit
/// cannot be audited after the fact.
pub struct Rng(pub [u64; 4]);

impl Rng {
    pub fn seed(s: u64) -> Rng {
        // splitmix64 to spread a single seed over the full state.
        let mut z = s;
        let mut next = || {
            z = z.wrapping_add(0x9E37_79B9_7F4A_7C15);
            let mut x = z;
            x = (x ^ (x >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
            x = (x ^ (x >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
            x ^ (x >> 31)
        };
        Rng([next(), next(), next(), next()])
    }

    pub fn next_u64(&mut self) -> u64 {
        let s = &mut self.0;
        let result = s[1].wrapping_mul(5).rotate_left(7).wrapping_mul(9);
        let t = s[1] << 17;
        s[2] ^= s[0];
        s[3] ^= s[1];
        s[1] ^= s[2];
        s[0] ^= s[3];
        s[2] ^= t;
        s[3] = s[3].rotate_left(45);
        result
    }

    /// Uniform in [0, 1).
    pub fn next_f64(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 * (1.0 / (1u64 << 53) as f64)
    }

    /// Uniform in [0, n).
    pub fn below(&mut self, n: u64) -> u64 {
        if n == 0 {
            return 0;
        }
        // Rejection sampling: taking a modulus directly would bias the low values.
        let zone = u64::MAX - u64::MAX % n;
        loop {
            let x = self.next_u64();
            if x < zone {
                return x % n;
            }
        }
    }
}

/// Build a matched case-control set.
///
/// `span` is the inclusive `(first, last)` day of the study period; controls are
/// only accepted inside it, so cases near either end simply contribute fewer
/// controls rather than being dropped or being given out-of-span partners.
///
/// Returns rows in case order, each case immediately followed by its controls.
pub fn build(
    quakes: &[Quake],
    grid: &Grid,
    scheme: Scheme,
    controls_per_case: usize,
    span: (f64, f64),
    rng: &mut Rng,
) -> Vec<Row> {
    let mut out = Vec::with_capacity(quakes.len() * (controls_per_case + 1));
    for (i, q) in quakes.iter().enumerate() {
        if q.day < span.0 || q.day > span.1 {
            continue;
        }
        let cell = grid.cell(q.lat_deg, q.lon_deg);
        let (lat, lon) = grid.centre(cell);
        let stratum = i as u32;
        out.push(Row {
            day: q.day,
            cell,
            lat_deg: lat,
            lon_deg: lon,
            case: true,
            stratum,
            magnitude: q.magnitude,
        });

        // The time-stratified referents are determined by the calendar rather than
        // drawn, so they are enumerated rather than sampled.
        if let Scheme::TimeStratified { spacing_days } = scheme {
            let step = spacing_days.max(1) as f64;
            let day_floor = q.day.floor();
            let frac = q.day - day_floor;
            let (y, m, _) = civil_from_days(day_floor as i64);
            let mut k = 1i64;
            // Walk outwards in both directions until the month is left behind.
            let mut refs: Vec<f64> = Vec::new();
            loop {
                let mut any = false;
                for sign in [-1.0, 1.0] {
                    let t = day_floor + sign * k as f64 * step;
                    let (ty, tm, _) = civil_from_days(t as i64);
                    if ty == y && tm == m {
                        any = true;
                        if t >= span.0 && t <= span.1 {
                            refs.push(t + frac);
                        }
                    }
                }
                if !any {
                    break;
                }
                k += 1;
            }
            for t in refs {
                out.push(Row {
                    day: t,
                    cell,
                    lat_deg: lat,
                    lon_deg: lon,
                    case: false,
                    stratum,
                    magnitude: f64::NAN,
                });
            }
            continue;
        }

        let mut used: Vec<f64> = Vec::with_capacity(controls_per_case);
        for _ in 0..controls_per_case {
            // A few attempts, then give up on this control rather than loop forever
            // when the span leaves no room -- a case on the first day of the study
            // has only forward offsets available.
            let mut placed = None;
            for _ in 0..32 {
                let t = match scheme {
                    Scheme::TimeStratified { .. } => unreachable!("handled above"),
                    Scheme::DayOffset { max_days } => {
                        let k = 1.0 + rng.below(max_days as u64) as f64;
                        let sign = if rng.next_u64() & 1 == 0 { -1.0 } else { 1.0 };
                        q.day + sign * k
                    }
                    Scheme::Window { window_days } => {
                        q.day + (rng.next_f64() * 2.0 - 1.0) * window_days
                    }
                    Scheme::Uniform => span.0 + rng.next_f64() * (span.1 - span.0),
                };
                if t < span.0 || t > span.1 {
                    continue;
                }
                if used.iter().any(|u| (u - t).abs() < 1e-9) {
                    continue;
                }
                placed = Some(t);
                break;
            }
            let Some(t) = placed else { continue };
            used.push(t);
            out.push(Row {
                day: t,
                cell,
                lat_deg: lat,
                lon_deg: lon,
                case: false,
                stratum,
                magnitude: f64::NAN,
            });
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn quake(day: f64, lat: f64, lon: f64) -> Quake {
        Quake { day, lat_deg: lat, lon_deg: lon, depth_km: 10.0, magnitude: 5.6 }
    }

    fn sample_quakes(n: usize, rng: &mut Rng) -> Vec<Quake> {
        (0..n)
            .map(|_| {
                quake(
                    rng.next_f64() * 17_800.0 - 8_766.0,
                    rng.next_f64() * 120.0 - 60.0,
                    rng.next_f64() * 360.0 - 180.0,
                )
            })
            .collect()
    }

    #[test]
    fn civil_dates_round_trip_against_known_values() {
        assert_eq!(civil_from_days(0), (2000, 1, 1));
        assert_eq!(civil_from_days(-1), (1999, 12, 31));
        assert_eq!(civil_from_days(59), (2000, 2, 29));   // 2000 was a leap year
        assert_eq!(civil_from_days(-8766), (1976, 1, 1));
        assert_eq!(civil_from_days(9131), (2024, 12, 31));
        assert_eq!(civil_from_days(9132), (2025, 1, 1));
    }

    #[test]
    fn time_stratified_referents_share_the_month_and_the_weekday() {
        let g = Grid::new(100.0);
        let mut rng = Rng::seed(21);
        let qs = sample_quakes(400, &mut rng);
        let span = (-8766.0, 9132.0);
        let rows = build(&qs, &g, Scheme::TimeStratified { spacing_days: 7 }, 0, span, &mut rng);

        let mut case_day = 0.0;
        let mut per_stratum = 0;
        let mut counts = Vec::new();
        for r in &rows {
            if r.case {
                if per_stratum > 0 {
                    counts.push(per_stratum);
                }
                per_stratum = 0;
                case_day = r.day;
                continue;
            }
            per_stratum += 1;
            let (cy, cm, _) = civil_from_days(case_day.floor() as i64);
            let (ry, rm, _) = civil_from_days(r.day.floor() as i64);
            assert_eq!((cy, cm), (ry, rm), "referent left the case's calendar month");
            let d = (r.day - case_day).abs();
            assert!((d / 7.0 - (d / 7.0).round()).abs() < 1e-9,
                    "offset {d} is not a multiple of 7 days");
            // Whole weeks, so both local solar time and day of week are held fixed.
            assert!((d - d.round()).abs() < 1e-9);
        }
        counts.push(per_stratum);
        assert!(counts.iter().all(|&c| (2..=4).contains(&c)),
                "referent counts out of range: {:?}",
                counts.iter().copied().collect::<std::collections::BTreeSet<_>>());
    }

    /// The property the whole scheme exists for.
    ///
    /// If the case is exchangeable with its referents then, under a uniform hazard,
    /// its rank among the stratum's times must be uniform. The bidirectional design
    /// fails this by construction: the case is always near the middle.
    #[test]
    fn the_case_is_exchangeable_with_its_referents() {
        let g = Grid::new(100.0);
        let span = (-8000.0, 9000.0);

        let rank_spread = |scheme: Scheme| -> f64 {
            let mut rng = Rng::seed(77);
            let qs = sample_quakes(6000, &mut rng);
            let rows = build(&qs, &g, scheme, 4, span, &mut rng);
            let mut mid = 0usize;
            let mut total = 0usize;
            let mut i = 0;
            while i < rows.len() {
                let case = rows[i];
                let mut times = vec![case.day];
                let mut j = i + 1;
                while j < rows.len() && !rows[j].case {
                    times.push(rows[j].day);
                    j += 1;
                }
                if times.len() >= 4 {
                    let n = times.len();
                    let rank = times.iter().filter(|&&t| t < case.day).count();
                    // Is the case in the interior rather than at either extreme?
                    if rank > 0 && rank < n - 1 {
                        mid += 1;
                    }
                    total += 1;
                }
                i = j;
            }
            mid as f64 / total as f64
        };

        // With n times and a uniformly-placed case, the chance of landing in the
        // interior is (n-2)/n: 0.5 at n = 4, 0.6 at n = 5.
        let ts = rank_spread(Scheme::TimeStratified { spacing_days: 7 });
        assert!(
            (0.4..=0.7).contains(&ts),
            "time-stratified case lands mid-stratum {ts:.3} of the time, expected ~0.5-0.6"
        );

        // The bidirectional design should visibly fail the same check.
        let bd = rank_spread(Scheme::DayOffset { max_days: 5 });
        assert!(
            bd > ts + 0.1,
            "expected the bidirectional design to centre the case ({bd:.3}) \
             far more than the time-stratified one ({ts:.3})"
        );
    }

    #[test]
    fn rng_is_reproducible_and_uniform() {
        let a: Vec<u64> = { let mut r = Rng::seed(42); (0..8).map(|_| r.next_u64()).collect() };
        let b: Vec<u64> = { let mut r = Rng::seed(42); (0..8).map(|_| r.next_u64()).collect() };
        assert_eq!(a, b, "same seed must give the same stream");
        let c: Vec<u64> = { let mut r = Rng::seed(43); (0..8).map(|_| r.next_u64()).collect() };
        assert_ne!(a, c);

        // Rough uniformity over ten buckets -- enough to catch a broken generator.
        let mut r = Rng::seed(7);
        let mut hist = [0usize; 10];
        let n = 200_000;
        for _ in 0..n {
            hist[(r.next_f64() * 10.0) as usize] += 1;
        }
        for (i, &h) in hist.iter().enumerate() {
            let dev = (h as f64 - n as f64 / 10.0).abs() / (n as f64 / 10.0);
            assert!(dev < 0.05, "bucket {i} off by {:.1}%", dev * 100.0);
        }
    }

    #[test]
    fn below_is_unbiased_near_the_modulus_boundary() {
        let mut r = Rng::seed(11);
        let mut hist = [0usize; 3];
        for _ in 0..120_000 {
            hist[r.below(3) as usize] += 1;
        }
        for (i, &h) in hist.iter().enumerate() {
            assert!((h as f64 / 40_000.0 - 1.0).abs() < 0.02, "value {i} count {h}");
        }
    }

    #[test]
    fn controls_share_their_cases_cell_and_never_its_time() {
        let g = Grid::new(100.0);
        let mut rng = Rng::seed(1);
        let qs = sample_quakes(500, &mut rng);
        let rows = build(&qs, &g, Scheme::DayOffset { max_days: 5 }, 10, (-8766.0, 9132.0), &mut rng);

        let mut i = 0;
        while i < rows.len() {
            let case = rows[i];
            assert!(case.case, "row {i} should start a stratum");
            let mut j = i + 1;
            while j < rows.len() && !rows[j].case {
                assert_eq!(rows[j].cell, case.cell, "control left the case's cell");
                assert_eq!(rows[j].stratum, case.stratum);
                assert!((rows[j].day - case.day).abs() >= 1.0 - 1e-9, "control at the case time");
                assert!(rows[j].magnitude.is_nan(), "control carries a magnitude");
                j += 1;
            }
            i = j;
        }
    }

    #[test]
    fn day_offsets_hold_local_solar_time_fixed() {
        // The defining property of the scheme: the fractional part of the day, which
        // is what sets local solar time, must be identical for case and controls.
        let g = Grid::new(100.0);
        let mut rng = Rng::seed(2);
        let qs = sample_quakes(300, &mut rng);
        let rows = build(&qs, &g, Scheme::DayOffset { max_days: 7 }, 6, (-8766.0, 9132.0), &mut rng);
        for w in rows.windows(2) {
            if !w[1].case {
                let d = (w[1].day - w[0].day).abs();
                assert!(
                    (d - d.round()).abs() < 1e-9,
                    "offset {d} is not a whole number of days"
                );
            }
        }
    }

    #[test]
    fn window_scheme_stays_inside_its_window() {
        let g = Grid::new(100.0);
        let mut rng = Rng::seed(3);
        let qs = sample_quakes(300, &mut rng);
        let span = (-8766.0, 9132.0);
        let rows = build(&qs, &g, Scheme::Window { window_days: 30.0 }, 8, span, &mut rng);
        let mut case_day = 0.0;
        for r in &rows {
            if r.case {
                case_day = r.day;
            } else {
                assert!((r.day - case_day).abs() <= 30.0 + 1e-9);
                assert!(r.day >= span.0 && r.day <= span.1);
            }
        }
    }

    #[test]
    fn cases_outside_the_span_are_excluded() {
        let g = Grid::new(100.0);
        let mut rng = Rng::seed(4);
        let qs = vec![quake(-20_000.0, 10.0, 10.0), quake(0.0, 10.0, 10.0), quake(20_000.0, 10.0, 10.0)];
        let rows = build(&qs, &g, Scheme::Uniform, 4, (-8766.0, 9132.0), &mut rng);
        let cases: Vec<&Row> = rows.iter().filter(|r| r.case).collect();
        assert_eq!(cases.len(), 1);
        assert_eq!(cases[0].day, 0.0);
    }

    #[test]
    fn a_case_at_the_edge_of_the_span_still_gets_controls() {
        // One-sided offsets must still work rather than the case silently losing
        // every control and quietly dropping out of the analysis.
        let g = Grid::new(100.0);
        let mut rng = Rng::seed(5);
        let span = (0.0, 1000.0);
        let qs = vec![quake(0.0, 5.0, 5.0), quake(1000.0, 5.0, 5.0)];
        let rows = build(&qs, &g, Scheme::DayOffset { max_days: 5 }, 6, span, &mut rng);
        for stratum in 0..2u32 {
            let n = rows.iter().filter(|r| r.stratum == stratum && !r.case).count();
            assert!(n >= 4, "edge case {stratum} got only {n} controls");
        }
        assert!(rows.iter().all(|r| r.day >= span.0 && r.day <= span.1));
    }

    /// The calibration test the whole design rests on.
    ///
    /// Events are generated with **no** dependence on anything astronomical — times
    /// uniform in the span, positions uniform on the sphere. Any statistic that
    /// separates cases from controls here is separating them by construction, which
    /// means the design leaks. This is the check that would have caught the
    /// block-shift and detection-bias errors recorded in `docs/07-research-log.md`
    /// before they produced a confident wrong answer.
    #[test]
    fn a_signal_free_catalogue_produces_no_case_control_separation() {
        let g = Grid::new(100.0);
        let span = (-8766.0, 9132.0);
        for scheme in [
            Scheme::DayOffset { max_days: 5 },
            Scheme::Window { window_days: 90.0 },
            Scheme::Uniform,
        ] {
            let mut rng = Rng::seed(99);
            let qs = sample_quakes(4000, &mut rng);
            let rows = build(&qs, &g, scheme, 8, span, &mut rng);

            // Test statistic: mean of cos and sin of the lunar-ish angle a naive
            // analysis would reach for -- the phase of a 29.53-day cycle. Under the
            // null the case mean and the control mean must agree to within the
            // sampling error of the smaller group.
            for (label, period) in [("synodic", 29.530_588), ("solar day", 1.0), ("annual", 365.25)] {
                let stat = |sel: bool| -> (f64, f64, usize) {
                    let v: Vec<f64> = rows
                        .iter()
                        .filter(|r| r.case == sel)
                        .map(|r| (std::f64::consts::TAU * r.day / period).cos())
                        .collect();
                    let n = v.len();
                    let m = v.iter().sum::<f64>() / n as f64;
                    let var = v.iter().map(|x| (x - m).powi(2)).sum::<f64>() / n as f64;
                    (m, var, n)
                };
                let (mc, vc, nc) = stat(true);
                let (mk, vk, nk) = stat(false);
                let se = (vc / nc as f64 + vk / nk as f64).sqrt();
                let z = (mc - mk) / se;
                assert!(
                    z.abs() < 4.0,
                    "{scheme:?} / {label}: cases and controls separate at z = {z:.2} \
                     under a signal-free catalogue"
                );
            }
        }
    }
}
