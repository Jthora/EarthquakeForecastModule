#!/usr/bin/env python3
"""Is the scan's top hit astronomy, or the calendar?

    python3 scripts/calendar_check.py --data /Volumes/.../m40_timestrat

Neptune-Pluto is the slowest-moving pair in the solar system: their relative
longitude drifts about 0.006 degrees a day. Any feature built from it is, to an
excellent approximation, a smooth monotonic function of the date. So if the scan
reports Neptune-Pluto aspects as its strongest hits, there are two readings, and
they are not close in plausibility:

  1. the Neptune-Pluto angle influences earthquakes
  2. the design has a residual dependence on WHEN, and the smoothest available
     function of when is picking it up

This distinguishes them without touching the feature matrix. Feature vectors are
not needed -- the day of every row is in the metadata, so the same conditional
score statistic can be computed for date itself, and for plain functions of it.
If z(date) is comparable to z(Neptune-Pluto), reading 2 is correct and the hit is
an artefact of the design rather than a finding about the sky.

Also reported: the same statistic under strata thinned to be mutually distant in
space and time. If the inflation is residual clustering -- nearby events sharing
nearly identical slow-feature values, which a within-stratum permutation cannot
see because it treats strata as independent -- thinning removes it.
"""

import argparse
import numpy as np

TEST_START = 6210.0


def score(values, starts, sizes, case_pos):
    """Conditional score z for one column, given stratum layout."""
    seg = np.repeat(np.arange(len(sizes)), sizes)
    m = np.add.reduceat(values, starts) / sizes
    c = values - m[seg]
    U = c[starts + case_pos].sum()
    sq = np.add.reduceat(c * c, starts) / sizes
    V = sq.sum()
    return U / np.sqrt(V) if V > 0 else 0.0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True)
    ap.add_argument("--perms", type=int, default=2000)
    ap.add_argument("--filter-col", default=None)
    ap.add_argument("--filter-val", default=None)
    a = ap.parse_args()

    rows = np.genfromtxt(a.data + ".rows.csv", delimiter=",", names=True,
                         dtype=None, encoding="utf-8")
    day = rows["day"].astype(np.float64)
    stratum = rows["stratum"]
    case = rows["case"].astype(bool)
    lat = rows["lat"].astype(np.float64)
    lon = rows["lon"].astype(np.float64)
    n = len(rows)

    starts = np.concatenate(([0], np.flatnonzero(np.diff(stratum)) + 1))
    sizes = np.diff(np.append(starts, n))
    case_pos = np.array([np.flatnonzero(case[s:s + z])[0] if case[s:s + z].sum() == 1 else -1
                         for s, z in zip(starts, sizes)])
    cd = np.array([day[s + p] if p >= 0 else np.nan for s, p in zip(starts, case_pos)])
    keep = (case_pos >= 0) & (sizes > 1) & (cd < TEST_START)
    if a.filter_col and a.filter_val:
        lab = np.array([str(rows[a.filter_col][s]) for s in starts])
        keep &= (lab == a.filter_val)
        print(f"restricted to {a.filter_col} == {a.filter_val}")
    print(f"{keep.sum()} usable strata")

    def subset(sel):
        idx = np.concatenate([np.arange(s, s + z) for s, z in zip(starts[sel], sizes[sel])])
        sz = sizes[sel]
        st = np.concatenate(([0], np.cumsum(sz)[:-1]))
        return idx, st, sz, case_pos[sel]

    idx, st, sz, cp = subset(keep)
    d = day[idx]

    print("\nconditional score z for plain functions of date "
          "(no feature matrix involved):")
    cols = {
        "date": d,
        "date^2": (d / 1000.0) ** 2,
        "Neptune-Pluto proxy: cos(2pi*date/181000)": np.cos(2 * np.pi * d / 181_000.0),
        "same, sin": np.sin(2 * np.pi * d / 181_000.0),
        "day of month": np.array([(x % 30.44) for x in d]),
    }
    for name, v in cols.items():
        print(f"  z = {score(v, st, sz, cp):+8.3f}   {name}")

    # Empirical null for these, by permuting which row of each stratum is the case.
    rng = np.random.default_rng(20260823)
    null = {k: [] for k in cols}
    for _ in range(a.perms):
        p = (rng.random(len(sz)) * sz).astype(np.int64)
        for k, v in cols.items():
            null[k].append(score(v, st, sz, p))
    print("\n  permutation null (same statistic, case reassigned at random):")
    for k in cols:
        arr = np.array(null[k])
        z = score(cols[k], st, sz, cp)
        pv = (np.sum(np.abs(arr) >= abs(z)) + 1) / (a.perms + 1)
        print(f"    {k:<44} null sd {arr.std():.3f}   p = {pv:.4f}")

    # How much thinning does it actually take? Discarding 91% of the catalogue is
    # a large price, so sweep the criterion and find the mildest one that removes
    # the artefact. z(date) is the readout: it should be consistent with N(0,1),
    # since the calendar cannot cause earthquakes.
    ks = np.flatnonzero(keep)
    order = np.argsort(cd[ks])
    print(f"\n{'km':>6} {'days':>6} {'strata':>8} {'z(date)':>9} {'z(N-P)':>9}"
          f" {'z(day of month)':>16}")
    for thin_km, thin_days in ((0, 0), (100, 30), (100, 180), (250, 180),
                               (250, 365), (500, 365), (1000, 365)):
        if thin_km == 0:
            sel = keep.copy()
        else:
            acc_lat, acc_lon, acc_day = [], [], []
            acc = []
            for si in ks[order]:
                t = cd[si]
                la, lo_ = lat[starts[si]], lon[starts[si]]
                ok = True
                for j in range(len(acc) - 1, -1, -1):
                    if t - acc_day[j] > thin_days:
                        break
                    dlat = (la - acc_lat[j]) * 111.19
                    dlon = (lo_ - acc_lon[j]) * 111.19 * np.cos(np.radians(la))
                    if dlat * dlat + dlon * dlon < thin_km ** 2:
                        ok = False
                        break
                if ok:
                    acc.append(si); acc_lat.append(la); acc_lon.append(lo_); acc_day.append(t)
            sel = np.zeros(len(sizes), dtype=bool)
            sel[np.array(acc)] = True
        idx2, st2, sz2, cp2 = subset(sel)
        d2 = day[idx2]
        zd = score(d2, st2, sz2, cp2)
        znp = score(np.cos(2 * np.pi * d2 / 181_000.0), st2, sz2, cp2)
        zdm = score(np.array([(x % 30.44) for x in d2]), st2, sz2, cp2)
        print(f"{thin_km:>6} {thin_days:>6} {int(sel.sum()):>8} "
              f"{zd:>+9.3f} {znp:>+9.3f} {zdm:>+16.3f}", flush=True)


if __name__ == "__main__":
    main()
