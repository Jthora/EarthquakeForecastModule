//! Equal-area spatial cells.
//!
//! Binning the globe by fixed degrees of latitude and longitude makes polar cells
//! vanishingly small — a 1°×1° cell at 80° N has a sixth the area of one at the
//! equator. That would show up in a model as a spurious latitude dependence, since
//! cell area sets the expected count directly.
//!
//! Here latitude bands have constant height and each band is divided into however
//! many longitude cells keep the area near the target:
//!
//! ```text
//! n_lon(band) = round( 2π (sin φ₂ − sin φ₁) / Δφ² )
//! ```
//!
//! which follows from the spherical zone area `2πR²(sin φ₂ − sin φ₁)` divided by the
//! target cell area `(RΔφ)²`. Bands near the poles collapse to a single cell rather
//! than to slivers.

/// Mean Earth radius, km.
pub const EARTH_RADIUS_KM: f64 = 6371.0;

/// An equal-area grid over the sphere.
#[derive(Debug, Clone)]
pub struct Grid {
    /// Band height in radians.
    dlat: f64,
    /// Longitude cells in each band, south to north.
    n_lon: Vec<usize>,
    /// Index of the first cell of each band; one longer than `n_lon`, so the last
    /// entry is the total cell count.
    offset: Vec<usize>,
}

impl Grid {
    /// A grid whose cells are about `target_km` on a side.
    pub fn new(target_km: f64) -> Grid {
        assert!(target_km > 0.0 && target_km < EARTH_RADIUS_KM, "implausible cell size");
        let dlat = target_km / EARTH_RADIUS_KM;
        let n_bands = (std::f64::consts::PI / dlat).round().max(1.0) as usize;
        let dlat = std::f64::consts::PI / n_bands as f64;

        let mut n_lon = Vec::with_capacity(n_bands);
        let mut offset = Vec::with_capacity(n_bands + 1);
        let mut total = 0usize;
        for i in 0..n_bands {
            let lat1 = -std::f64::consts::FRAC_PI_2 + i as f64 * dlat;
            let lat2 = lat1 + dlat;
            let n = (std::f64::consts::TAU * (lat2.sin() - lat1.sin()) / (dlat * dlat))
                .round()
                .max(1.0) as usize;
            offset.push(total);
            total += n;
            n_lon.push(n);
        }
        offset.push(total);
        Grid { dlat, n_lon, offset }
    }

    /// Total number of cells.
    pub fn len(&self) -> usize {
        self.offset[self.offset.len() - 1]
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Number of latitude bands.
    pub fn bands(&self) -> usize {
        self.n_lon.len()
    }

    /// Cell containing a point. Latitude is clamped to the poles; longitude wraps.
    pub fn cell(&self, lat_deg: f64, lon_deg: f64) -> usize {
        let lat = lat_deg.to_radians().clamp(-std::f64::consts::FRAC_PI_2, std::f64::consts::FRAC_PI_2);
        let band = (((lat + std::f64::consts::FRAC_PI_2) / self.dlat).floor() as isize)
            .clamp(0, self.n_lon.len() as isize - 1) as usize;
        let n = self.n_lon[band];
        let mut frac = (lon_deg / 360.0).fract();
        if frac < 0.0 {
            frac += 1.0;
        }
        let col = ((frac * n as f64).floor() as usize).min(n - 1);
        self.offset[band] + col
    }

    /// Centre of a cell, in degrees.
    pub fn centre(&self, cell: usize) -> (f64, f64) {
        let band = self.band_of(cell);
        let lat = -std::f64::consts::FRAC_PI_2 + (band as f64 + 0.5) * self.dlat;
        let col = cell - self.offset[band];
        let lon = (col as f64 + 0.5) / self.n_lon[band] as f64 * 360.0;
        (lat.to_degrees(), if lon > 180.0 { lon - 360.0 } else { lon })
    }

    /// Area of a cell, km².
    pub fn area_km2(&self, cell: usize) -> f64 {
        let band = self.band_of(cell);
        let lat1 = -std::f64::consts::FRAC_PI_2 + band as f64 * self.dlat;
        let lat2 = lat1 + self.dlat;
        std::f64::consts::TAU * EARTH_RADIUS_KM * EARTH_RADIUS_KM * (lat2.sin() - lat1.sin())
            / self.n_lon[band] as f64
    }

    fn band_of(&self, cell: usize) -> usize {
        assert!(cell < self.len(), "cell {cell} out of range");
        // offset is sorted; partition_point gives the first band starting past cell.
        self.offset.partition_point(|&o| o <= cell) - 1
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn every_point_lands_in_a_valid_cell() {
        let g = Grid::new(100.0);
        for lat in [-90.0, -89.999, -45.0, 0.0, 33.3, 89.999, 90.0] {
            for lon in [-180.0, -179.9, -77.0, 0.0, 139.7, 179.9, 180.0, 360.0, -360.0] {
                let c = g.cell(lat, lon);
                assert!(c < g.len(), "lat {lat} lon {lon} gave cell {c} of {}", g.len());
            }
        }
    }

    #[test]
    fn cell_areas_are_nearly_equal_everywhere() {
        // The whole point of the scheme. A fixed-degree grid would fail this badly.
        let g = Grid::new(100.0);
        let areas: Vec<f64> = (0..g.len()).map(|c| g.area_km2(c)).collect();
        let target = 100.0 * 100.0;
        let (min, max) = areas.iter().fold((f64::MAX, 0.0f64), |(lo, hi), &a| (lo.min(a), hi.max(a)));
        assert!(
            max / min < 1.6,
            "areas span {min:.0} to {max:.0} km2, ratio {:.2}",
            max / min
        );
        let mean = areas.iter().sum::<f64>() / areas.len() as f64;
        assert!((mean / target - 1.0).abs() < 0.1, "mean area {mean:.0} vs target {target:.0}");
    }

    #[test]
    fn areas_sum_to_the_surface_of_the_sphere() {
        // Catches gaps and overlaps in one shot: if bands double-counted or missed
        // any zone, the total would not close.
        let g = Grid::new(250.0);
        let total: f64 = (0..g.len()).map(|c| g.area_km2(c)).sum();
        let sphere = 4.0 * std::f64::consts::PI * EARTH_RADIUS_KM * EARTH_RADIUS_KM;
        assert!(
            (total / sphere - 1.0).abs() < 1e-12,
            "cells cover {total:.0} km2, sphere is {sphere:.0} km2"
        );
    }

    #[test]
    fn centres_map_back_to_their_own_cells() {
        let g = Grid::new(150.0);
        for c in (0..g.len()).step_by(7) {
            let (lat, lon) = g.centre(c);
            assert_eq!(g.cell(lat, lon), c, "cell {c} centre ({lat}, {lon}) round-trip");
        }
    }

    #[test]
    fn polar_bands_collapse_to_three_cells() {
        // The polar zone's area tends to pi*R^2*dlat^2/... precisely pi times the
        // target cell area as the band height shrinks, so the cap always divides
        // into three wedges regardless of resolution. Worth pinning: it is the
        // scheme refusing to make slivers, and it is resolution-independent.
        for km in [50.0, 100.0, 250.0] {
            let g = Grid::new(km);
            assert_eq!(g.n_lon[0], 3, "south cap at {km} km");
            assert_eq!(g.n_lon[g.bands() - 1], 3, "north cap at {km} km");
        }
        let g = Grid::new(100.0);
        // And the equatorial band should have the full circumference's worth.
        let mid = g.bands() / 2;
        let expect = (std::f64::consts::TAU * EARTH_RADIUS_KM / 100.0).round() as usize;
        assert!(
            (g.n_lon[mid] as isize - expect as isize).abs() <= 2,
            "equatorial band has {} cells, expected about {expect}",
            g.n_lon[mid]
        );
    }

    #[test]
    fn adjacent_longitudes_land_in_adjacent_cells_and_wrap() {
        let g = Grid::new(100.0);
        let lat = 0.2;
        let a = g.cell(lat, 179.99);
        let b = g.cell(lat, -179.99);
        // Crossing the date line must wrap within the band, not jump bands.
        let band_a = g.offset.partition_point(|&o| o <= a) - 1;
        let band_b = g.offset.partition_point(|&o| o <= b) - 1;
        assert_eq!(band_a, band_b, "date line crossing changed band");
        assert_ne!(a, b);
    }
}
