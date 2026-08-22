#!/usr/bin/env python3
"""Conditional logistic regression on the matched design.

    python3 scripts/train_conditional_logit.py --data /Volumes/.../m55_dayoffset

Model A of docs/22-model-grid.md. Fits on 1976-2009, selects the L2 strength on
2010-2016, and does NOT touch 2017-2024 -- opening the test set is a separate
script run once, deliberately, after this one has named a winner.

The likelihood is conditional on the matched set: within a stratum of one case
and k controls the model says which row is the case, so it is a softmax over the
stratum and the intercept cancels. That is what makes the design immune to the
case rate -- no calibration of absolute probability is attempted or needed.

Standardisation is folded into the coefficients rather than applied to the data.
With gamma = beta/sigma,

    eta_i = sum_k (x_ik - mu_k)/sigma_k * beta_k = (X @ gamma)_i - mu . gamma

and the mu term is constant within a stratum, so it cancels in the softmax. The
fit therefore runs on the raw matrix with no standardised copy and no per-chunk
arithmetic; sigma survives only in the penalty, as lambda/2 * sum (gamma_k sigma_k)^2.
"""

import argparse, json, os, sys, time
import numpy as np
from scipy.optimize import minimize

CHUNK = 4096


def load(prefix):
    names = [l.strip() for l in open(prefix + ".names") if l.strip()]
    rows = np.genfromtxt(prefix + ".rows.csv", delimiter=",", names=True,
                         dtype=None, encoding="utf-8")
    n, d = len(rows), len(names)
    size = os.path.getsize(prefix + ".f32")
    expect = n * d * 4
    if size != expect:
        sys.exit(f"matrix is {size} bytes, expected {expect} for {n} x {d}")
    X = np.memmap(prefix + ".f32", dtype=np.float32, mode="r", shape=(n, d))
    return names, rows, X


def segments(stratum):
    """Start index of each stratum, given rows grouped by stratum."""
    change = np.flatnonzero(np.diff(stratum)) + 1
    return np.concatenate(([0], change))


def eta_of(X, gamma, idx):
    """X[idx] @ gamma, in chunks so no large temporary is ever built."""
    out = np.empty(len(idx), dtype=np.float64)
    g32 = gamma.astype(np.float32)
    for a in range(0, len(idx), CHUNK):
        b = min(a + CHUNK, len(idx))
        out[a:b] = X[idx[a:b]] @ g32
    return out


def xt_dot(X, v, idx, d):
    """X[idx].T @ v, accumulated in float64."""
    acc = np.zeros(d, dtype=np.float64)
    for a in range(0, len(idx), CHUNK):
        b = min(a + CHUNK, len(idx))
        acc += X[idx[a:b]].T @ v[a:b].astype(np.float32)
    return acc


class Split:
    """Rows of one period, grouped into strata."""

    def __init__(self, X, rows, mask, d):
        idx = np.flatnonzero(mask)
        order = np.argsort(rows["stratum"][idx], kind="stable")
        self.idx = idx[order]
        self.stratum = rows["stratum"][self.idx]
        self.y = rows["case"][self.idx].astype(np.float64)
        self.starts = segments(self.stratum)
        self.seg = np.repeat(np.arange(len(self.starts)),
                             np.diff(np.append(self.starts, len(self.idx))))
        self.sizes = np.diff(np.append(self.starts, len(self.idx)))
        self.X, self.d = X, d
        # A stratum whose case fell outside the period, or which lost every
        # control, carries no information and would divide by zero.
        keep = np.isin(self.seg, np.flatnonzero((self.sizes > 1)))
        keep &= np.isin(self.seg, np.flatnonzero(
            np.bincount(self.seg, weights=self.y, minlength=len(self.sizes)) == 1))
        if not keep.all():
            self.idx, self.y = self.idx[keep], self.y[keep]
            self.stratum = self.stratum[keep]
            self.starts = segments(self.stratum)
            self.seg = np.repeat(np.arange(len(self.starts)),
                                 np.diff(np.append(self.starts, len(self.idx))))
            self.sizes = np.diff(np.append(self.starts, len(self.idx)))

    @property
    def n_strata(self):
        return len(self.sizes)

    def logp_case(self, gamma):
        """Log probability assigned to the true case in each stratum."""
        eta = eta_of(self.X, gamma, self.idx)
        m = np.maximum.reduceat(eta, self.starts)
        e = np.exp(eta - m[self.seg])
        denom = np.add.reduceat(e, self.starts)
        eta_case = eta[self.y == 1]
        return eta_case - m - np.log(denom), eta, e, denom

    def neg_loglik_and_grad(self, gamma, lam, sigma2, active):
        lp, eta, e, denom = self.logp_case(gamma)
        nll = -lp.sum()
        p = e / denom[self.seg]
        g = xt_dot(self.X, p - self.y, self.idx, self.d)
        pen = lam * gamma * sigma2
        obj = nll + 0.5 * lam * float(np.sum(sigma2 * gamma * gamma))
        grad = g + pen
        grad[~active] = 0.0
        return obj, grad

    def info_gain(self, gamma):
        """Bits per event above a model that assigns 1/size to every row."""
        lp, _, _, _ = self.logp_case(gamma)
        null = -np.log(self.sizes.astype(np.float64))
        return float(np.mean(lp - null) / np.log(2))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True)
    ap.add_argument("--out", default=None)
    ap.add_argument("--maxiter", type=int, default=500)
    ap.add_argument("--lambdas", default="1e-4,1e-3,1e-2,1e-1,1,10,100,1000")
    a = ap.parse_args()
    out = a.out or (a.data + ".logit.json")

    names, rows, X = load(a.data)
    d = len(names)
    print(f"{len(rows)} rows x {d} features")

    # Split by the date of each stratum's CASE, so a case and its controls never
    # straddle a boundary -- controls can sit up to a year away under the window
    # scheme, and splitting on the row's own date would leak across periods.
    case_day = {}
    for s, day, is_case in zip(rows["stratum"], rows["day"], rows["case"]):
        if is_case:
            case_day[s] = day
    day_of_stratum = np.array([case_day.get(s, np.nan) for s in rows["stratum"]])

    TRAIN_END, VAL_END = 3652.0, 6210.0     # 2010-01-01 and 2017-01-01
    tr = day_of_stratum < TRAIN_END
    va = (day_of_stratum >= TRAIN_END) & (day_of_stratum < VAL_END)
    te = day_of_stratum >= VAL_END
    train, val = Split(X, rows, tr, d), Split(X, rows, va, d)
    print(f"train {train.n_strata} strata, validate {val.n_strata} strata, "
          f"test {int(np.sum(te & (rows['case'] == 1)))} strata (sealed)")

    print("computing training-period feature scales...")
    t0 = time.time()
    n = len(train.idx)
    s1 = np.zeros(d); s2 = np.zeros(d)
    for i in range(0, n, CHUNK):
        blk = X[train.idx[i:min(i + CHUNK, n)]].astype(np.float64)
        s1 += blk.sum(0); s2 += (blk * blk).sum(0)
    mu = s1 / n
    var = np.maximum(s2 / n - mu * mu, 0.0)
    sigma = np.sqrt(var)
    active = sigma > 1e-9
    print(f"  {int((~active).sum())} constant columns dropped, "
          f"{int(active.sum())} active ({time.time()-t0:.0f}s)")
    sigma2 = np.where(active, var, 1.0)

    results = []
    for lam in [float(x) for x in a.lambdas.split(",")]:
        t0 = time.time()
        g0 = np.zeros(d)
        r = minimize(
            lambda g: train.neg_loglik_and_grad(g, lam, sigma2, active),
            g0, jac=True, method="L-BFGS-B",
            options={"maxiter": a.maxiter, "gtol": 1e-6, "maxcor": 10},
        )
        ig_tr, ig_va = train.info_gain(r.x), val.info_gain(r.x)
        results.append({"lambda": lam, "train_ig": ig_tr, "val_ig": ig_va,
                        "iters": int(r.nit), "converged": bool(r.success),
                        "gnorm": float(np.max(np.abs(r.jac)))})
        print(f"  lambda {lam:<8g}  train {ig_tr:+.5f}  validate {ig_va:+.5f} bits"
              f"   {r.nit} iters, {time.time()-t0:.0f}s")
        np.save(a.data + f".gamma_lam{lam:g}.npy", r.x)

    best = max(results, key=lambda r: r["val_ig"])
    print(f"\nbest on validation: lambda {best['lambda']:g}, "
          f"{best['val_ig']:+.5f} bits/event")
    med = float(np.median([r["val_ig"] for r in results]))
    print(f"median across {len(results)} configs: {med:+.5f} bits/event")
    json.dump({"data": a.data, "results": results, "best": best, "median": med,
               "n_train": train.n_strata, "n_val": val.n_strata,
               "active_features": int(active.sum())},
              open(out, "w"), indent=2)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
