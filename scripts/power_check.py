#!/usr/bin/env python3
"""How large an effect would this design actually detect?

    python3 scripts/power_check.py --data ~/eqf-work/m55_dayoffset

A null result means nothing without this. With 9781 features and 8414 training
strata, a pipeline that cannot recover a *planted* effect of realistic size is
indistinguishable from one where no effect exists -- and the distinction is the
whole point.

So: keep the real feature matrix exactly as it is, throw away the real labels,
and generate new ones from a known effect. Within each stratum the case is drawn
with probability proportional to exp(beta . x), which is exactly the model being
fitted. Then sweep beta and see where the validation information gain rises out
of the noise.

Two shapes of planted signal, because they are not equally easy to find:

  concentrated   all of the effect on a single feature
  distributed    the same total effect spread over 20 features

Real signal, if it existed, would more likely be distributed -- an aspect that
matters would matter across several harmonics and frames at once.

Effect sizes are quoted as the odds multiplier per standard deviation of the
feature, so beta = 0.05 is a 5% modulation of relative rate per SD. Published
tidal-triggering effects on ordinary crust sit at a few percent, and this
programme previously bounded them below 3.88% at M2 (docs/19-results.md).
"""

import argparse, os, sys, time
import numpy as np
from scipy.optimize import minimize

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import importlib.util
spec = importlib.util.spec_from_file_location(
    "cl", os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "train_conditional_logit.py"))
cl = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cl)


def plant(split, beta_vec, rng):
    """Redraw which row of each stratum is the case, from a known coefficient."""
    eta = (split.Xm @ beta_vec.astype(np.float32)).astype(np.float64)
    y = np.zeros(len(eta))
    for s in range(len(split.starts)):
        lo = split.starts[s]
        hi = lo + split.sizes[s]
        w = np.exp(eta[lo:hi] - eta[lo:hi].max())
        y[lo + rng.choice(split.sizes[s], p=w / w.sum())] = 1.0
    return y


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True)
    ap.add_argument("--max-controls", type=int, default=4)
    ap.add_argument("--lam", type=float, default=1000.0,
                    help="the L2 strength that won on validation for the real fit")
    ap.add_argument("--maxiter", type=int, default=200)
    ap.add_argument("--betas", default="0,0.02,0.05,0.1,0.2,0.5")
    a = ap.parse_args()

    names, rows, X = cl.load(a.data)
    d = len(names)
    case_day = {s: dd for s, dd, c in zip(rows["stratum"], rows["day"], rows["case"]) if c}
    dos = np.array([case_day.get(s, np.nan) for s in rows["stratum"]])
    tr = dos < 3652.0
    va = (dos >= 3652.0) & (dos < 6210.0)
    train = cl.Split(X, rows, tr, d, a.max_controls)
    val = cl.Split(X, rows, va, d, a.max_controls)
    train.load(); val.load()

    n = len(train.idx)
    s1 = np.zeros(d); s2 = np.zeros(d)
    for i in range(0, n, cl.CHUNK):
        blk = train.Xm[i:min(i + cl.CHUNK, n)].astype(np.float64)
        s1 += blk.sum(0); s2 += (blk * blk).sum(0)
    mu = s1 / n
    var = np.maximum(s2 / n - mu * mu, 0.0)
    sigma = np.sqrt(var)
    active = sigma > 1e-9
    mu32 = mu.astype(np.float32)
    inv32 = np.where(active, 1.0 / np.where(active, sigma, 1.0), 0.0).astype(np.float32)
    train.standardise(mu32, inv32)
    val.standardise(mu32, inv32)
    print(f"{train.n_strata} train strata, {val.n_strata} validate, "
          f"{int(active.sum())} active features, lambda {a.lam:g}")

    # Features to carry the planted signal: a lunar one for the concentrated
    # case, and a spread of aspect harmonics for the distributed case. Chosen by
    # name so the check is reproducible and not cherry-picked after the fact.
    idx_of = {nm: i for i, nm in enumerate(names)}
    conc = idx_of.get("geo.moon.syn.h2.cos")
    dist = [i for nm, i in idx_of.items()
            if nm.startswith("geo.asp.") and nm.endswith(".h2.cos")][:20]
    print(f"concentrated on '{names[conc]}'; distributed over {len(dist)} h2 aspects\n")

    rng = np.random.default_rng(20260822)
    print(f"{'beta':>6}  {'shape':<12} {'train IG':>10} {'valid IG':>10}   verdict")
    for b in [float(x) for x in a.betas.split(",")]:
        for shape in ("concentrated", "distributed"):
            bv = np.zeros(d)
            if b > 0:
                if shape == "concentrated":
                    bv[conc] = b
                else:
                    # Same total signal-to-noise: spreading over k features with
                    # each coefficient b/sqrt(k) keeps ||beta|| fixed.
                    for i in dist:
                        bv[i] = b / np.sqrt(len(dist))
            y_tr = plant(train, bv, rng)
            y_va = plant(val, bv, rng)
            old_tr, old_va = train.y, val.y
            train.y, val.y = y_tr, y_va
            r = minimize(lambda g: train.neg_loglik_and_grad(g, a.lam, active),
                         np.zeros(d), jac=True, method="L-BFGS-B",
                         options={"maxiter": a.maxiter, "gtol": 1e-8, "ftol": 1e-14})
            ig_tr, ig_va = train.info_gain(r.x), val.info_gain(r.x)
            train.y, val.y = old_tr, old_va
            verdict = "DETECTED" if ig_va >= 0.01 else ("weak" if ig_va > 0 else "invisible")
            print(f"{b:>6.2f}  {shape:<12} {ig_tr:>+10.5f} {ig_va:>+10.5f}   {verdict}")
            if b == 0:
                break   # no shape distinction at zero


if __name__ == "__main__":
    main()
