#!/usr/bin/env python3
"""Null calibration for the training pipeline (docs/21-preregistration.md section 9).

    python3 scripts/test_pipeline_null.py

Three checks, all of which must pass before any real fit is trusted:

  1. no signal      synthetic features independent of the label -> validation IG
                    must not be reliably positive
  2. known signal   one feature genuinely predictive -> must be recovered, so a
                    null result cannot be blamed on a model that cannot learn
  3. label shuffle  real structure, labels permuted within stratum -> IG ~ 0

Check 2 matters as much as the other two. A pipeline that reports zero because
it is broken is indistinguishable, from the outside, from one that reports zero
because there is nothing there.
"""

import os, sys, tempfile
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import importlib.util
spec = importlib.util.spec_from_file_location(
    "cl", os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "train_conditional_logit.py"))
cl = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cl)

RNG = np.random.default_rng(20260822)
K = 10          # controls per case
D = 60          # features
N_STRATA = 3000


def make(prefix, signal_beta=0.0, shuffle=False):
    """Write a synthetic dataset in the featurise output format."""
    n = N_STRATA * (K + 1)
    X = RNG.standard_normal((n, D)).astype(np.float32)
    stratum = np.repeat(np.arange(N_STRATA), K + 1)
    case = np.zeros(n, dtype=np.int64)

    if signal_beta == 0.0:
        # The case is a uniformly random member of its stratum: nothing in X
        # carries any information about which row it is.
        for s in range(N_STRATA):
            case[s * (K + 1) + RNG.integers(K + 1)] = 1
    else:
        # The case is drawn with probability proportional to exp(beta * x[:, 0]),
        # which is exactly the conditional logistic model the fitter assumes.
        for s in range(N_STRATA):
            lo = s * (K + 1)
            w = np.exp(signal_beta * X[lo:lo + K + 1, 0].astype(np.float64))
            case[lo + RNG.choice(K + 1, p=w / w.sum())] = 1

    if shuffle:
        for s in range(N_STRATA):
            lo = s * (K + 1)
            blk = case[lo:lo + K + 1].copy()
            RNG.shuffle(blk)
            case[lo:lo + K + 1] = blk

    # Spread strata across the split boundaries used by the trainer.
    day = np.repeat(RNG.uniform(-8766.0, 9132.0, N_STRATA), K + 1)

    X.tofile(prefix + ".f32")
    with open(prefix + ".names", "w") as f:
        for i in range(D):
            f.write(f"f{i}\n")
    with open(prefix + ".rows.csv", "w") as f:
        f.write("day,cell,lat,lon,case,stratum,magnitude\n")
        for i in range(n):
            f.write(f"{day[i]:.9f},0,0.0,0.0,{case[i]},{stratum[i]},"
                    f"{'5.5' if case[i] else ''}\n")
    return prefix


def fit(prefix, lam):
    names, rows, X = cl.load(prefix)
    d = len(names)
    case_day = {s: dd for s, dd, c in zip(rows["stratum"], rows["day"], rows["case"]) if c}
    dos = np.array([case_day.get(s, np.nan) for s in rows["stratum"]])
    tr, va = dos < 3652.0, (dos >= 3652.0) & (dos < 6210.0)
    train, val = cl.Split(X, rows, tr, d), cl.Split(X, rows, va, d)
    n = len(train.idx)
    blk = X[train.idx].astype(np.float64)
    var = np.maximum(blk.var(0), 0.0)
    active = np.sqrt(var) > 1e-9
    from scipy.optimize import minimize
    r = minimize(lambda g: train.neg_loglik_and_grad(g, lam, np.where(active, var, 1.0), active),
                 np.zeros(d), jac=True, method="L-BFGS-B",
                 options={"maxiter": 500, "gtol": 1e-8})
    return train.info_gain(r.x), val.info_gain(r.x), r.x, train, val


def main():
    tmp = tempfile.mkdtemp(prefix="nullcal-")
    fail = []

    print("1. no signal -- validation IG must not be reliably positive")
    vals = []
    for rep in range(5):
        p = make(os.path.join(tmp, f"null{rep}"), signal_beta=0.0)
        ig_tr, ig_va, _, _, _ = fit(p, lam=1.0)
        vals.append(ig_va)
        print(f"   rep {rep}: train {ig_tr:+.5f}  validate {ig_va:+.5f} bits")
    mean = float(np.mean(vals))
    print(f"   mean validation IG {mean:+.5f} bits")
    if mean > 0.01:
        fail.append(f"null data gave {mean:+.5f} bits/event, above the 0.01 success bar")
    # Overfitting on train while validation stays flat is the expected signature
    # and is not itself a failure -- it is what the held-out period is for.

    print("\n2. known signal -- must be recovered")
    p = make(os.path.join(tmp, "sig"), signal_beta=0.5)
    ig_tr, ig_va, gamma, _, _ = fit(p, lam=1.0)
    print(f"   train {ig_tr:+.5f}  validate {ig_va:+.5f} bits")
    print(f"   coefficient on the signal feature: {gamma[0]:+.4f} "
          f"(others |max| {np.abs(gamma[1:]).max():.4f})")
    if ig_va < 0.02:
        fail.append(f"planted signal recovered at only {ig_va:+.5f} bits -- fitter is broken")
    if abs(gamma[0]) < 3 * np.abs(gamma[1:]).max():
        fail.append("signal feature does not stand out from the noise features")

    print("\n3. label shuffle within stratum -- IG must collapse")
    p = make(os.path.join(tmp, "shuf"), signal_beta=0.5, shuffle=True)
    ig_tr, ig_va, _, _, _ = fit(p, lam=1.0)
    print(f"   train {ig_tr:+.5f}  validate {ig_va:+.5f} bits")
    if ig_va > 0.01:
        fail.append(f"shuffled labels still gave {ig_va:+.5f} bits/event")

    print()
    if fail:
        for f in fail:
            print(f"FAIL: {f}")
        sys.exit(1)
    print("null calibration passed: the pipeline finds a planted signal and "
          "finds nothing in noise.")


if __name__ == "__main__":
    main()
