# 20 — Data Acquisition Plan

Everything needed to feed the model. Sizes below are verified, not estimated.

**Guiding rule:** no mechanistic filtering. If a quantity can be computed or
downloaded, it goes in. The discipline lives in validation, not feature selection.

---

## Tier 0 — already in hand

| Data | Size | Status |
|---|---|---|
| DE440s ephemeris + kernels | 33 MB | ✅ |
| ComCat M4.0+ 1976–2024 | 488,214 events | ✅ |
| ComCat M5.5+ 1970–2024 | 25,962 events | ✅ |
| GCMT focal mechanisms | 67,263 solutions | ✅ |
| Parkfield LFEs | 1,528,117 | ✅ |
| Cascadia tremor | 678,084 | ✅ |
| Apollo deep moonquakes | 6,954 | ✅ |
| SPOTL + GOT4.7/FES2004 ocean models | 213 MB | ✅ |
| Ocean loading at 18,316 event sites | M2, O1 | ✅ |

---

## Tier 1 — computed, no download

Everything ephemeris-derived. This is a build task, not an acquisition task.

**Three frames** — geocentric, heliocentric, barycentric — for 13 bodies (Sun,
Moon, Mercury…Pluto, mean and true node):

- ecliptic longitude, latitude, distance
- **signed** speed (retrograde is negative), and acceleration
- declination and right ascension
- tropical **and** sidereal longitude (they differ by precession, ~24° and drifting)
- galactic coordinates and galactic-centre angle

**Derived per epoch:**

| Family | Content |
|---|---|
| Aspects | every pair, every frame, Fourier-encoded `cos nΔλ, sin nΔλ` for n = 1…24 |
| Declination aspects | parallels and contraparallels |
| Chart shape | circular variance, largest gap, concentration, cluster count, centroid — the classical bundle/bowl/splash shapes as continuous numbers |
| Whole-chart resonance | CosmicCypher's base-N aggregate, **with base 7 and 11 fixed** and cosine restored |
| Eclipses | proximity, magnitude, type — syzygy and node alignment together |
| Lunar | distance, declination, libration, node position, phase angle |
| Stations | proximity to retrograde/direct turning, via `ph_core::events` |
| Tidal | full constituent set, ΔCFS, dΔCFS/dt, tensor eigenvalues, principal-axis direction |
| Commensurabilities | the d'Alembert multi-body set, filtered to resolvable periods |
| Site-local | local apparent sidereal time, ASC/MC, solar altitude (day/night) |

### Storage: keep primitives, derive features

Naive materialisation is **~12 GB**: 430k hourly epochs × ~3,500 features × 8 bytes.

Don't. Store the **~200 primitives** per epoch — positions, speeds, declinations
across three frames — which is **690 MB**, and derive aspects, harmonics and shapes
in the training loop. They are cheap trigonometry on the primitives.

Twenty-fold saving, and the derived set can change without recomputing ephemerides.

---

## Tier 2 — free, direct, no credentials

All verified reachable and sized.

| Source | Content | Cadence | Span | Size |
|---|---|---|---|---|
| **OMNI2** (NASA SPDF) | solar wind speed, density, temperature; IMF Bx/By/Bz; plus Kp, Dst, AE | **hourly** | 1963– | **184 MB** |
| **GFZ Potsdam** | Kp, ap, Ap, sunspot number, F10.7 — one file | 3-hourly | 1932– | 5.5 MB |
| **IERS** `finals2000A.all.csv` | **LOD, polar motion x/y, UT1−UTC** | daily | 1973– | 4.0 MB |
| **SILSO** (Royal Obs. Belgium) | daily sunspot number v2.0 | daily | 1818– | 2.9 MB |

```
https://spdf.gsfc.nasa.gov/pub/data/omni/low_res_omni/omni2_all_years.dat
https://kp.gfz-potsdam.de/app/files/Kp_ap_Ap_SN_F107_since_1932.txt
https://datacenter.iers.org/data/csv/finals2000A.all.csv
https://www.sidc.be/SILSO/DATA/SN_d_tot_V2.0.csv
```

**OMNI2 is the prize** — one hourly file carrying solar wind, interplanetary
magnetic field and the geomagnetic indices together, back to 1963. It alone covers
most of the space-weather brainstorm.

Total Tier 2: **~200 MB**, four downloads, no accounts.

---

## Tier 3 — free but credential-gated

Registration only, no cost. Worth doing; not blocking.

| Source | Content | Gate |
|---|---|---|
| GRACE / GRACE-FO mascons (JPL) | mass loading, monthly, 2002– | NASA Earthdata login |
| GLDAS (NASA) | land surface hydrology, 3-hourly | NASA Earthdata login |
| ERA5 (Copernicus) | surface pressure, 3 TB if global — subset to cells | CDS API key |

**Note the span limit:** GRACE starts in 2002, so mass-loading features cover only
~40% of the catalogue. Either restrict those runs or let the model handle missing
values explicitly.

---

## Tier 4 — patchy, best effort

Schumann resonance has no well-maintained long archive. Tomsk and HeartMath series
exist with gaps. Try, accept partial coverage, and do not block on it.

---

## Time alignment

Everything resamples onto a **common hourly grid**, 1976–2024 ≈ 430,000 rows.

| Native cadence | Handling |
|---|---|
| ephemeris | evaluated directly at grid times |
| hourly (OMNI2) | direct |
| 3-hourly (Kp) | step interpolation — the index *is* a 3-hour bin |
| daily (IERS, SILSO) | linear interpolation |
| monthly (GRACE) | linear, and flagged as interpolated |

⚠ Interpolating slow series creates artificial smoothness. Carry an explicit
`is_interpolated` flag per family so the model can learn to distrust it, rather
than us pretending the resolution is real.

---

## Acquisition gotchas already paid for

Both of these silently corrupted a catalogue during earlier work:

1. **PNSN caps responses at 20,000 events** with HTTP 200 and no truncation flag. A
   yearly request returned exactly 20,000 and looked fine.
2. **ComCat's cap bites below M5.0**, so chunk yearly rather than by decade.
3. **Empty windows return 404, not an empty result**, so `curl --fail` aborts a run
   on legitimately-empty early years.

Every fetch script asserts expected row counts and warns near any cap.

---

## Order

```
1. Tier 1 primitives          the build; nothing blocks it
2. Tier 2 downloads           ~200 MB, four files, parallel with 1
3. Feature derivation layer   primitives -> ~3,500 features
4. Time alignment             one hourly table
5. Tier 3 (registration)      GRACE/GLDAS/ERA5 as a second pass
```

Tier 1 is the long pole and needs no network. Tier 2 can run alongside it in
minutes.
