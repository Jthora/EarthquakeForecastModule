#!/usr/bin/env python3
"""Conditional score test, one feature at a time, with a permutation null.

    python3 scripts/score_scan.py --data ~/eqf-work/m55_dayoffset

The penalised 9781-feature fit turned out to be powerless: `power_check.py`
showed that even a planted effect of 0.5 log-odds per SD -- a 50% modulation,
far larger than anything in the tidal-triggering literature -- is invisible to
it. With more features than strata, L2 must shrink so hard that real signal goes
with the noise. A null result from that model is a statement about the model.

This is the powerful test. Under the matched null the case is equally likely to
be any row of its stratum, so for feature k

    U_k = sum_s ( x[case] - mean_s(x) )
    V_k = sum_s ( mean_s(x^2) - mean_s(x)^2 )
    z_k = U_k / sqrt(V_k)

is standard normal, exactly, with no model fitted and nothing to converge. It
costs two passes over the data and has the full power of the design behind each
feature individually.

The price is multiplicity: 9781 tests. The threshold is set by permutation
rather than by Bonferroni, because the features are heavily correlated -- cos
and sin of neighbouring harmonics of the same pair are nearly the same
question -- so the effective number of independent tests is far below 9781 and
Bonferroni would be needlessly conservative. Permuting which row of each stratum
is the case preserves the correlation structure exactly and gives the true null
distribution of max |z|.

The test period is not touched: the scan runs on 1976-2016 only.
"""

import argparse, json, os, sys, time
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import importlib.util
spec = importlib.util.spec_from_file_location(
    "cl", os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "train_conditional_logit.py"))
cl = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cl)

TEST_START = 6210.0     # 2017-01-01


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True)
    ap.add_argument("--max-controls", type=int, default=4)
    ap.add_argument("--perms", type=int, default=500)
    ap.add_argument("--top", type=int, default=25)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    out = a.out or (a.data + ".scan.json")

    names, rows, X = cl.load(a.data)
    d = len(names)
    case_day = {s: dd for s, dd, c in zip(rows["stratum"], rows["day"], rows["case"]) if c}
    dos = np.array([case_day.get(s, np.nan) for s in rows["stratum"]])
    pre = dos < TEST_START
    sp = cl.Split(X, rows, pre, d, a.max_controls)
    t0 = time.time()
    gb = sp.load()
    print(f"{sp.n_strata} strata, {len(sp.idx)} rows, {d} features "
          f"({gb:.2f} GB, {time.time()-t0:.0f}s)")
    print(f"period: day {dos[pre].min():.0f} to {dos[pre].max():.0f} (test sealed)")

    # Centre every feature within its stratum, in place. After this the null has
    # mean zero by construction and the case rows are all that matter.
    t0 = time.time()
    n_s = sp.sizes
    sums = np.add.reduceat(sp.Xm, sp.starts, axis=0)
    means = (sums / n_s[:, None]).astype(np.float32)
    sq = np.add.reduceat(sp.Xm.astype(np.float32) ** 2, sp.starts, axis=0)
    var_s = np.maximum(sq / n_s[:, None] - means.astype(np.float64) ** 2, 0.0)
    V = var_s.sum(0)
    for i in range(len(sp.starts)):
        lo, hi = sp.starts[i], sp.starts[i] + n_s[i]
        sp.Xm[lo:hi] -= means[i]
    print(f"centred within strata ({time.time()-t0:.0f}s)")

    case_rows = np.flatnonzero(sp.y == 1)
    U = sp.Xm[case_rows].sum(0, dtype=np.float64)
    ok = V > 1e-12
    z = np.zeros(d)
    z[ok] = U[ok] / np.sqrt(V[ok])
    print(f"{int(ok.sum())} testable features")

    # Is the observed z field even standard normal? If it is wider than N(0,1),
    # something in the design is correlating rows within a stratum -- residual
    # aftershock clustering being the obvious candidate -- and every p-value
    # computed from it would be inflated.
    zo = z[ok]
    print(f"\nz distribution: mean {zo.mean():+.4f}, sd {zo.std():.4f} "
          f"(expected 0.000, 1.000)")
    print(f"  |z| > 2: {int((np.abs(zo) > 2).sum())} of {int(ok.sum())} "
          f"({100*np.mean(np.abs(zo) > 2):.2f}%, expect 4.55%)")
    print(f"  |z| > 3: {int((np.abs(zo) > 3).sum())} "
          f"({100*np.mean(np.abs(zo) > 3):.3f}%, expect 0.27%)")
    print(f"  max |z| = {np.abs(zo).max():.3f}")

    # Permutation null for max |z|, preserving the correlation between features.
    t0 = time.time()
    rng = np.random.default_rng(20260822)
    offs = np.repeat(sp.starts, 1)
    maxes = np.empty(a.perms)
    for p in range(a.perms):
        pick = offs + (rng.random(len(n_s)) * n_s).astype(np.int64)
        Up = sp.Xm[pick].sum(0, dtype=np.float64)
        zp = np.zeros(d)
        zp[ok] = Up[ok] / np.sqrt(V[ok])
        maxes[p] = np.abs(zp[ok]).max()
        if p % 100 == 0:
            print(f"  permutation {p}/{a.perms} ({time.time()-t0:.0f}s)", flush=True)
    thresh = np.quantile(maxes, 0.95)
    obs = np.abs(zo).max()
    pval = float((np.sum(maxes >= obs) + 1) / (a.perms + 1))
    print(f"\npermutation null for max |z| over {a.perms} draws:")
    print(f"  median {np.median(maxes):.3f}, 95th percentile {thresh:.3f}")
    print(f"  observed {obs:.3f}  ->  family-wise p = {pval:.4f}")
    # How many independent tests does that imply? Useful as a sanity number.
    from scipy.stats import norm
    eff = np.log(0.5) / np.log(2 * norm.cdf(-np.median(maxes)) / 2 + 1e-300)
    print(f"  implies roughly {eff:.0f} effectively independent tests "
          f"(of {int(ok.sum())} features)")

    order = np.argsort(-np.abs(z))
    print(f"\ntop {a.top} features by |z|:")
    for i in order[:a.top]:
        mark = " *" if abs(z[i]) >= thresh else ""
        print(f"  {z[i]:+8.3f}  {names[i]}{mark}")

    json.dump({"data": a.data, "n_strata": int(sp.n_strata),
               "n_features": int(ok.sum()),
               "z_mean": float(zo.mean()), "z_sd": float(zo.std()),
               "max_abs_z": float(obs), "perm_threshold_95": float(thresh),
               "family_wise_p": pval, "perms": a.perms,
               "top": [{"name": names[i], "z": float(z[i])} for i in order[:200]]},
              open(out, "w"), indent=2)
    print(f"\nwrote {out}")
    if pval < 0.05:
        print("At least one feature separates cases from controls beyond chance.")
    else:
        print("No feature separates cases from controls beyond chance.")


if __name__ == "__main__":
    main()
