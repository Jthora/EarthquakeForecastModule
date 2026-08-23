//! Build the design matrix: every matched row joined to its full feature vector.
//!
//!     cargo run --release -p eqf-dataset --bin featurise -- \
//!         --min-mag 5.5 --out /Volumes/2TB_EXT_1B/eqf-data/m55
//!
//! Writes three files next to `--out`:
//!
//! ```text
//! <out>.f32      row-major matrix, n_rows x n_features, little-endian f32
//! <out>.names    one feature name per line, in column order
//! <out>.rows.csv day, cell, lat, lon, case, stratum, magnitude
//! ```
//!
//! Raw f32 rather than a framed format because the consumer is a training loop
//! that memory-maps it and never parses anything. f32 rather than f64 because the
//! features are trigonometric quantities in [-1, 1] whose seventh digit is far
//! below anything a model can use, and the matrix is large enough that halving it
//! decides whether it fits in memory.
//!
//! # Batching
//!
//! Ephemeris lookups are an order of magnitude cheaper in batches than one at a
//! time, so rows are processed in chunks and each chunk makes one SPICE call per
//! frame rather than one per row.

use eqf::{comcat, gcmt};
use eqf_dataset::{cells::Grid, decluster, sampling, strata};
use ph_core::{chart, chart_cycles, chart_features, chart_local};
use rustspice_core::KernelSet;
use std::io::{BufWriter, Write};

const CHUNK: usize = 4096;

struct Config {
    min_mag: f64,
    max_harmonic: usize,
    max_base: usize,
    local_harmonic: usize,
    controls: usize,
    max_days: u32,
    scheme: String,
    window_days: f64,
    spacing_days: u32,
    block_years: i64,
    seed: u64,
    out: String,
    catalogue: String,
    limit: usize,
    format: String,
}

fn parse_args() -> Config {
    let mut c = Config {
        min_mag: 5.5,
        max_harmonic: 24,
        max_base: 24,
        local_harmonic: 12,
        controls: 10,
        max_days: 5,
        scheme: "dayoffset".into(),
        window_days: 90.0,
        spacing_days: 7,
        block_years: 6,
        seed: 20260822,
        out: "dataset".into(),
        catalogue: "data/comcat/global_m40.csv".into(),
        limit: usize::MAX,
        format: "comcat".into(),
    };
    let args: Vec<String> = std::env::args().skip(1).collect();
    let mut i = 0;
    // A closure would have to borrow `i` mutably while the match still reads it,
    // so the advance is written out.
    macro_rules! val {
        () => {{
            i += 1;
            args.get(i)
                .cloned()
                .unwrap_or_else(|| panic!("missing value for {}", args[i - 1]))
        }};
    }
    while i < args.len() {
        match args[i].as_str() {
            "--min-mag" => c.min_mag = val!().parse().unwrap(),
            "--max-harmonic" => c.max_harmonic = val!().parse().unwrap(),
            "--max-base" => c.max_base = val!().parse().unwrap(),
            "--local-harmonic" => c.local_harmonic = val!().parse().unwrap(),
            "--controls" => c.controls = val!().parse().unwrap(),
            "--max-days" => c.max_days = val!().parse().unwrap(),
            "--scheme" => c.scheme = val!(),
            "--window-days" => c.window_days = val!().parse().unwrap(),
            "--spacing-days" => c.spacing_days = val!().parse().unwrap(),
            "--block-years" => c.block_years = val!().parse().unwrap(),
            "--seed" => c.seed = val!().parse().unwrap(),
            "--out" => c.out = val!(),
            "--catalogue" => c.catalogue = val!(),
            "--format" => c.format = val!(),
            "--limit" => c.limit = val!().parse().unwrap(),
            other => panic!("unknown argument {other}"),
        }
        i += 1;
    }
    c
}

/// Station timing for an arbitrary instant, from a precomputed station list.
///
/// [`chart_cycles::station_timing`] works on a chart series; event times are not on
/// any grid, so stations are located once on a daily grid over the whole span and
/// each row is then placed against that list by binary search.
fn timing_at(stations: &[Vec<f64>], day: f64, span: (f64, f64)) -> Vec<(String, f64)> {
    let mut out = Vec::new();
    for (bi, times) in stations.iter().enumerate() {
        let body = chart::BODIES[bi];
        let k = times.partition_point(|&t| t <= day);
        let since = if k > 0 { day - times[k - 1] } else { day - span.0 };
        let until = if k < times.len() { times[k] - day } else { span.1 - day };
        let span_len = since + until;
        out.push((format!("geo.stn.since.{body}"), since));
        out.push((format!("geo.stn.until.{body}"), until));
        out.push((format!("geo.stn.nearest.{body}"), since.min(until)));
        out.push((
            format!("geo.stn.frac.{body}"),
            if span_len > 0.0 { since / span_len } else { 0.0 },
        ));
    }
    out
}

fn main() -> rustspice_core::Result<()> {
    let cfg = parse_args();
    let t_start = std::time::Instant::now();

    // GCMT carries the focal mechanisms, so when it is the source it is also the
    // event list -- associating two catalogues event by event would introduce
    // matching errors for no gain.
    let text = std::fs::read_to_string(&cfg.catalogue).expect("read catalogue");
    let (all, labels): (Vec<comcat::Quake>, Option<Vec<(&'static str, &'static str)>>) =
        if cfg.format == "gcmt" {
            let cmts = gcmt::parse_ndk(&text);
            let mut qs = Vec::with_capacity(cmts.len());
            let mut ls = Vec::with_capacity(cmts.len());
            for c in &cmts {
                qs.push(comcat::Quake {
                    day: c.day,
                    lat_deg: c.lat_deg,
                    lon_deg: c.lon_deg,
                    depth_km: c.depth_km,
                    magnitude: c.mw,
                });
                ls.push((
                    strata::mechanism(c).name(),
                    if strata::depth_class(c) == strata::Depth::Shallow { "shallow" } else { "deep" },
                ));
            }
            (qs, Some(ls))
        } else {
            (comcat::parse_catalog(&text), None)
        };
    // Labels must follow their events through declustering, so they are carried on
    // a key rather than an index -- gardner_knopoff returns a filtered copy.
    let label_of: std::collections::HashMap<u64, (&'static str, &'static str)> = match &labels {
        Some(ls) => all
            .iter()
            .zip(ls)
            .map(|(q, l)| (q.day.to_bits(), *l))
            .collect(),
        None => Default::default(),
    };
    let (span_lo, span_hi) = all
        .iter()
        .fold((f64::MAX, f64::MIN), |(a, b), q| (a.min(q.day), b.max(q.day)));
    let span = (span_lo, span_hi);

    let sel: Vec<_> = all.iter().copied().filter(|q| q.magnitude >= cfg.min_mag).collect();
    let (main, dep) = decluster::gardner_knopoff(&sel);
    println!(
        "catalogue {}: {} events at M{}+, {} independent ({} dependent removed)",
        cfg.catalogue,
        sel.len(),
        cfg.min_mag,
        main.len(),
        dep.len()
    );

    let grid = Grid::new(100.0);
    let scheme = match cfg.scheme.as_str() {
        "dayoffset" => sampling::Scheme::DayOffset { max_days: cfg.max_days },
        "window" => sampling::Scheme::Window { window_days: cfg.window_days },
        "uniform" => sampling::Scheme::Uniform,
        "timestratified" => sampling::Scheme::TimeStratified { spacing_days: cfg.spacing_days },
        "yearstratified" => sampling::Scheme::YearStratified { block_years: cfg.block_years },
        other => panic!("unknown scheme {other}"),
    };
    let mut rng = sampling::Rng::seed(cfg.seed);
    let mut rows = sampling::build(&main, &grid, scheme, cfg.controls, span, &mut rng);
    if rows.len() > cfg.limit {
        rows.truncate(cfg.limit);
    }
    let n_cases = rows.iter().filter(|r| r.case).count();
    println!(
        "design: {scheme:?}, {} controls per case -> {} rows ({n_cases} cases)",
        cfg.controls,
        rows.len()
    );

    let mut ks = KernelSet::new();
    for k in ["naif0012.tls", "de440s.bsp", "pck00011.tpc"] {
        ks.add_file(format!("kernels/{k}"))?;
    }
    let mut spice = ks.open()?;
    let epoch2000 = spice.parse_time("2000-01-01T00:00:00")?;

    // Stations, once, on a daily grid across the whole span.
    println!("locating stations on a daily grid...");
    let daily: Vec<f64> = {
        let n = (span.1 - span.0).ceil() as usize + 1;
        (0..n).map(|i| span.0 + i as f64).collect()
    };
    let daily_charts = chart::charts(&mut spice, &daily, chart::Frame::Geocentric, epoch2000)?;
    let stations_flat = chart_cycles::stations(&daily_charts);
    let mut stations: Vec<Vec<f64>> = vec![Vec::new(); chart::BODIES.len()];
    for s in &stations_flat {
        stations[s.body].push(s.day);
    }
    drop(daily_charts);
    println!("  {} stations across {} bodies", stations_flat.len(), chart::BODIES.len());

    let frames = [chart::Frame::Geocentric, chart::Frame::Heliocentric, chart::Frame::Barycentric];
    let mut matrix = BufWriter::with_capacity(1 << 22, std::fs::File::create(format!("{}.f32", cfg.out))?);
    let mut names_written: Option<Vec<String>> = None;
    let mut n_features = 0usize;

    let mut meta = BufWriter::new(std::fs::File::create(format!("{}.rows.csv", cfg.out))?);
    writeln!(meta, "day,cell,lat,lon,case,stratum,magnitude,mech,depthclass").unwrap();
    // A stratum's label is its case's label; controls inherit it so that filtering
    // on the column keeps whole strata rather than slicing them.
    let case_label: std::collections::HashMap<u32, (&str, &str)> = rows
        .iter()
        .filter(|r| r.case)
        .map(|r| {
            (
                r.stratum,
                *label_of.get(&r.day.to_bits()).unwrap_or(&("", "")),
            )
        })
        .collect();
    for r in &rows {
        let (mech, dep) = case_label.get(&r.stratum).copied().unwrap_or(("", ""));
        writeln!(
            meta,
            "{:.9},{},{:.4},{:.4},{},{},{},{},{}",
            r.day,
            r.cell,
            r.lat_deg,
            r.lon_deg,
            if r.case { 1 } else { 0 },
            r.stratum,
            if r.magnitude.is_nan() { String::new() } else { format!("{:.2}", r.magnitude) },
            mech,
            dep
        )
        .unwrap();
    }
    meta.flush().unwrap();

    let mut done = 0usize;
    for chunk in rows.chunks(CHUNK) {
        let days: Vec<f64> = chunk.iter().map(|r| r.day).collect();
        let mut per_frame = Vec::with_capacity(frames.len());
        for f in frames {
            per_frame.push(chart::charts(&mut spice, &days, f, epoch2000)?);
        }

        for (i, row) in chunk.iter().enumerate() {
            let mut fs = chart_features::FeatureSet::default();
            for (fi, frame_charts) in per_frame.iter().enumerate() {
                let _ = fi;
                fs.extend(chart_features::all(&frame_charts[i], cfg.max_harmonic, cfg.max_base));
            }
            let cyc = chart_cycles::all(&per_frame[0][i], cfg.max_harmonic);
            for (n, v) in cyc.names.iter().zip(&cyc.values) {
                fs.push(format!("geo.{n}"), *v);
            }
            for (n, v) in timing_at(&stations, row.day, span) {
                fs.push(n, v);
            }
            let site = chart_local::Site::from_degrees(row.lat_deg, row.lon_deg);
            fs.extend(chart_local::all(&per_frame[0][i], site, cfg.local_harmonic));

            if names_written.is_none() {
                n_features = fs.len();
                let mut f = std::fs::File::create(format!("{}.names", cfg.out))?;
                for n in &fs.names {
                    writeln!(f, "{n}").unwrap();
                }
                println!("{n_features} features per row");
                names_written = Some(fs.names.clone());
            }
            assert_eq!(fs.len(), n_features, "feature count changed at row {done}");

            // Non-finite values would poison training silently; fail loudly instead,
            // naming the feature, so the cause is fixed rather than imputed over.
            for (n, v) in fs.names.iter().zip(&fs.values) {
                assert!(v.is_finite(), "row {done}: {n} = {v}");
            }
            let bytes: Vec<u8> = fs.values.iter().flat_map(|v| (*v as f32).to_le_bytes()).collect();
            matrix.write_all(&bytes).unwrap();
            done += 1;
        }
        let frac = done as f64 / rows.len() as f64;
        let elapsed = t_start.elapsed().as_secs_f64();
        eprint!(
            "\r  {done}/{} rows  {:.1}%  {:.0}s elapsed, {:.0}s left   ",
            rows.len(),
            frac * 100.0,
            elapsed,
            elapsed / frac.max(1e-9) - elapsed
        );
    }
    matrix.flush().unwrap();
    eprintln!();

    let bytes = done as u64 * n_features as u64 * 4;
    println!(
        "wrote {}.f32: {done} rows x {n_features} features = {:.2} GB in {:.0}s",
        cfg.out,
        bytes as f64 / 1e9,
        t_start.elapsed().as_secs_f64()
    );
    Ok(())
}
