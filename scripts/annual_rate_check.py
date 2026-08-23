#!/usr/bin/env python3
"""Is global seismicity Poisson from year to year?

    python3 scripts/annual_rate_check.py

The year-stratified design compares an earthquake against the same calendar date
in other years of its block. That is only valid if the event was equally likely
to have fallen in any of those years -- if the annual rate is constant. If it is
not, every slow-moving feature is confounded with the rate history, because a
slow feature is essentially a smooth function of the year.

The design showed z(date) = -6.19, and thinning to 1000 km and 365 days barely
moved it (-3.82). Local clustering cannot explain that: thinning removes it. A
global year-to-year variation in rate would, because it moves every cell at once
and no amount of spatial separation touches it.

So the question is whether declustered annual counts are Poisson. If they are
over-dispersed, the year-stratified design is not repairable by any referent
choice, and slow configurations are untestable against this catalogue rather
than merely unsupported by it.
"""

import numpy as np
from scipy import stats
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import subprocess

def load_counts(min_mag):
    """Declustered annual counts, via the Rust inspector's own declustering."""
    out = subprocess.run(
        ["./target/release/annual_counts", str(min_mag)],
        capture_output=True, text=True)
    if out.returncode != 0:
        sys.exit(out.stderr[-2000:])
    years, counts = [], []
    for line in out.stdout.strip().splitlines():
        y, c = line.split()
        years.append(int(y)); counts.append(int(c))
    return np.array(years), np.array(counts)


def main():
    for mag in (4.0, 5.0, 5.5, 6.0):
        years, counts = load_counts(mag)
        # Use only fully-covered years.
        keep = (years >= 1976) & (years <= 2024)
        years, counts = years[keep], counts[keep]
        m = counts.mean()
        v = counts.var(ddof=1)
        # Poisson dispersion test: sum (x - m)^2 / m ~ chi2 with n-1 df.
        chi2 = ((counts - m) ** 2).sum() / m
        df = len(counts) - 1
        p = stats.chi2.sf(chi2, df)
        print(f"M{mag}+  {len(counts)} years, mean {m:8.1f}/yr, "
              f"variance/mean = {v/m:6.2f}   chi2 = {chi2:8.1f} on {df} df, "
              f"p = {p:.3e}")
        if p < 0.001:
            # How much of a slow feature's apparent signal could this explain?
            # A feature that is monotonic in year correlates with the rate history
            # at exactly the strength of that history's own variation.
            r = np.corrcoef(years, counts)[0, 1]
            print(f"        over-dispersed by {v/m:.1f}x Poisson; "
                  f"linear trend with year r = {r:+.3f}")
    print()
    print("Over-dispersion here means the year-stratified referent design cannot")
    print("be repaired: no choice of referent controls a rate that varies across")
    print("the very years being compared. Slow configurations are untestable")
    print("against this catalogue, which is a weaker statement than unsupported.")


if __name__ == "__main__":
    main()
