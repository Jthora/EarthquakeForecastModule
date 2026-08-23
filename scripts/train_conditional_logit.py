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

Features are standardised in place after loading, using training-period
statistics only.

An earlier version folded the standardisation into the coefficients instead:
with gamma = beta/sigma the mean term is constant within a stratum and cancels
in the softmax, so the fit could run on the raw matrix with sigma surviving only
in the penalty. That is algebraically correct and numerically useless. It leaves
L-BFGS looking at raw columns whose variances span 0.5 for a cosine to 4e8 for a
lunar distance in km, a condition number around 1e9, and the quasi-Newton
approximation cannot cope: at lambda = 1e-4, where a model with 9781 features
and 8414 strata should overfit to near log2(5) bits, it reached 0.0005 and
stopped. The synthetic null missed it because those features were all N(0, 1).
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


def eta_chunked(X, beta, idx, mu, sigma):
    """X[idx] standardised @ beta, straight from the memmap.

    For splits not held in RAM. The same training-period mu and sigma are applied
    here as were applied in place to the training split.
    """
    out = np.empty(len(idx), dtype=np.float64)
    b32 = beta.astype(np.float32)
    for a in range(0, len(idx), CHUNK):
        b = min(a + CHUNK, len(idx))
        blk = (X[idx[a:b]] - mu) * sigma
        out[a:b] = blk @ b32
    return out


def materialise(X, idx):
    """Copy the split's rows out of the memmap into a contiguous array.

    Without this the fit is disk-bound: L-BFGS touches every row twice per
    iteration, and re-reading 3.6 GB through a memmap costs ~37 s per iteration
    at 14% CPU. Held in RAM the same iteration is under a second. The cost is
    that the split must fit -- 92k rows x 9816 f32 is 3.6 GB -- so this is done
    per split rather than for the whole matrix.
    """
    out = np.empty((len(idx), X.shape[1]), dtype=np.float32)
    for a in range(0, len(idx), CHUNK):
        b = min(a + CHUNK, len(idx))
        out[a:b] = X[idx[a:b]]
    return out


class Split:
    """Rows of one period, grouped into strata."""

    def __init__(self, X, rows, mask, d, max_controls=None):
        idx = np.flatnonzero(mask)
        order = np.argsort(rows["stratum"][idx], kind="stable")
        idx = idx[order]
        if max_controls is not None:
            # Keep the case and the first `max_controls` controls of each stratum.
            # They were drawn in generator order, so taking a prefix is an unbiased
            # subsample, and it is deterministic given the dataset's seed.
            keep = np.zeros(len(idx), dtype=bool)
            seen = {}
            for j, i in enumerate(idx):
                s_id = rows["stratum"][i]
                if rows["case"][i]:
                    keep[j] = True
                else:
                    c = seen.get(s_id, 0)
                    if c < max_controls:
                        keep[j] = True
                        seen[s_id] = c + 1
            idx = idx[keep]
        self.idx = idx
        self.stratum = rows["stratum"][self.idx]
        self.y = rows["case"][self.idx].astype(np.float64)
        self.starts = segments(self.stratum)
        self.seg = np.repeat(np.arange(len(self.starts)),
                             np.diff(np.append(self.starts, len(self.idx))))
        self.sizes = np.diff(np.append(self.starts, len(self.idx)))
        self.X, self.d = X, d
        self.Xm = None
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

    def load(self):
        """Bring this split into RAM. Returns GB used."""
        if self.Xm is None:
            self.Xm = materialise(self.X, self.idx)
        return self.Xm.nbytes / 1e9

    def standardise(self, mu, inv_sigma):
        """Apply (x - mu) / sigma in place, in chunks, so no copy is made."""
        for a in range(0, len(self.Xm), CHUNK):
            b = min(a + CHUNK, len(self.Xm))
            self.Xm[a:b] -= mu
            self.Xm[a:b] *= inv_sigma
        self.standardised = True

    @property
    def n_strata(self):
        return len(self.sizes)

    def logp_case(self, beta, mu=None, inv_sigma=None):
        """Log probability assigned to the true case in each stratum."""
        if self.Xm is not None:
            eta = (self.Xm @ beta.astype(np.float32)).astype(np.float64)
        else:
            eta = eta_chunked(self.X, beta, self.idx, mu, inv_sigma)
        m = np.maximum.reduceat(eta, self.starts)
        e = np.exp(eta - m[self.seg])
        denom = np.add.reduceat(e, self.starts)
        eta_case = eta[self.y == 1]
        return eta_case - m - np.log(denom), eta, e, denom

    def neg_loglik_and_grad(self, beta, lam, active):
        lp, eta, e, denom = self.logp_case(beta)
        nll = -lp.sum()
        p = e / denom[self.seg]
        g = (self.Xm.T @ (p - self.y).astype(np.float32)).astype(np.float64)
        obj = nll + 0.5 * lam * float(beta @ beta)
        grad = g + lam * beta
        grad[~active] = 0.0
        return obj, grad

    def info_gain(self, beta, mu=None, inv_sigma=None):
        """Bits per event above a model that assigns 1/size to every row."""
        lp, _, _, _ = self.logp_case(beta, mu, inv_sigma)
        null = -np.log(self.sizes.astype(np.float64))
        return float(np.mean(lp - null) / np.log(2))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True)
    ap.add_argument("--out", default=None)
    ap.add_argument("--maxiter", type=int, default=500)
    ap.add_argument("--lambdas", default="1e-4,1e-3,1e-2,1e-1,1,10,100,1000")
    ap.add_argument("--max-controls", type=int, default=None,
                    help="keep at most this many controls per stratum (memory)")
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
    train = Split(X, rows, tr, d, a.max_controls)
    val = Split(X, rows, va, d, a.max_controls)
    print(f"train {train.n_strata} strata, validate {val.n_strata} strata, "
          f"test {int(np.sum(te & (rows['case'] == 1)))} strata (sealed)")
    print(f"controls per stratum: {float(np.mean(train.sizes)) - 1:.2f} (train)")

    # Only the training split goes into RAM: L-BFGS touches it twice per
    # iteration, whereas validation is scored once per lambda and can stream.
    t0 = time.time()
    gb = train.load()
    print(f"loaded train into RAM: {gb:.2f} GB ({time.time()-t0:.0f}s); "
          f"validate streams from disk")

    print("computing training-period feature scales...")
    t0 = time.time()
    n = len(train.idx)
    s1 = np.zeros(d); s2 = np.zeros(d)
    for i in range(0, n, CHUNK):
        blk = train.Xm[i:min(i + CHUNK, n)].astype(np.float64)
        s1 += blk.sum(0); s2 += (blk * blk).sum(0)
    mu = s1 / n
    var = np.maximum(s2 / n - mu * mu, 0.0)
    sigma = np.sqrt(var)
    active = sigma > 1e-9
    print(f"  {int((~active).sum())} constant columns dropped, "
          f"{int(active.sum())} active ({time.time()-t0:.0f}s)")
    # A constant column gets inv_sigma 0, so it standardises to exactly zero and
    # contributes nothing rather than dividing by zero.
    mu32 = mu.astype(np.float32)
    inv_sigma32 = np.where(active, 1.0 / np.where(active, sigma, 1.0), 0.0).astype(np.float32)
    t0 = time.time()
    train.standardise(mu32, inv_sigma32)
    print(f"  standardised training split in place ({time.time()-t0:.0f}s)")
    # Sanity: standardised columns must have unit variance, or the fit is being
    # handed the same ill-conditioned problem under a different name.
    chk = train.Xm[:8192, active].astype(np.float64)
    sd = chk.std(0)
    print(f"  column sd after standardising: min {sd.min():.3f}, "
          f"median {np.median(sd):.3f}, max {sd.max():.3f}")

    results = []
    for lam in [float(x) for x in a.lambdas.split(",")]:
        t0 = time.time()
        g0 = np.zeros(d)
        r = minimize(
            lambda g: train.neg_loglik_and_grad(g, lam, active),
            g0, jac=True, method="L-BFGS-B",
            options={"maxiter": a.maxiter, "gtol": 1e-8, "ftol": 1e-14,
                     "maxcor": 20},
        )
        ig_tr = train.info_gain(r.x)
        ig_va = val.info_gain(r.x, mu32, inv_sigma32)
        results.append({"lambda": lam, "train_ig": ig_tr, "val_ig": ig_va,
                        "iters": int(r.nit), "converged": bool(r.success),
                        "gnorm": float(np.max(np.abs(r.jac)))})
        print(f"  lambda {lam:<8g}  train {ig_tr:+.5f}  validate {ig_va:+.5f} bits"
              f"   {r.nit} iters, {time.time()-t0:.0f}s")
        np.save(a.data + f".beta_lam{lam:g}.npy", r.x)
    np.save(a.data + ".scale.npy", np.vstack([mu, sigma]))

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
