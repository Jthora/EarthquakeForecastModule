#!/usr/bin/env python3
"""Model B: gradient-boosted trees on the matched design.

    ~/eqf-work/venv/bin/python scripts/train_gbt.py --data ~/eqf-work/m55_dayoffset

Model A (penalised conditional logistic) is linear and, at 9781 features against
8414 strata, provably powerless: power_check.py could not recover a planted
effect of 0.5 log-odds per SD. Trees are the reason to expect better. They pick
a few features per split instead of fitting every coefficient at once, and they
can represent interactions -- an aspect that only matters when the Moon is also
near a node -- which no additive model can express.

The objective is the same conditional likelihood, so the information gains are
directly comparable between the two models. Within a stratum of one case and k
controls the scores go through a softmax and the model is asked which row is the
case, giving

    grad_i = p_i - y_i          hess_i = p_i (1 - p_i)

Trees are scale-invariant, so no standardisation is applied here -- the raw
matrix is binned by LightGBM directly.
"""

import argparse, gc, json, os, sys, time
import numpy as np
import lightgbm as lgb

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import importlib.util
spec = importlib.util.spec_from_file_location(
    "cl", os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "train_conditional_logit.py"))
cl = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cl)

TRAIN_END, VAL_END = 3652.0, 6210.0


def softmax_by_stratum(preds, starts, sizes, seg):
    m = np.maximum.reduceat(preds, starts)
    e = np.exp(preds - m[seg])
    denom = np.add.reduceat(e, starts)
    return e / denom[seg]


def make_objective(starts, sizes, seg, y):
    def obj(preds, _dset):
        p = softmax_by_stratum(preds, starts, sizes, seg)
        grad = p - y
        hess = np.maximum(p * (1.0 - p), 1e-6)
        return grad, hess
    return obj


def make_eval(splits):
    """One evaluator for every dataset.

    LightGBM applies each feval to *all* valid_sets, so an evaluator that closed
    over a single split's stratum boundaries would be handed the other split's
    predictions. Dispatching on row count keeps one function correct for both;
    the two splits differ in length, which is asserted at construction.
    """
    by_len = {}
    for starts, sizes, seg, y in splits:
        n = int(np.sum(sizes))
        assert n not in by_len, "splits have equal row counts; cannot dispatch"
        by_len[n] = (starts, sizes, seg, y, -np.log(sizes.astype(np.float64)))

    def ev(preds, _dset):
        starts, sizes, seg, y, null = by_len[len(preds)]
        m = np.maximum.reduceat(preds, starts)
        e = np.exp(preds - m[seg])
        denom = np.add.reduceat(e, starts)
        lp = preds[y == 1] - m - np.log(denom)
        return "ig", float(np.mean(lp - null) / np.log(2)), True
    return ev


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True)
    ap.add_argument("--max-controls", type=int, default=4)
    ap.add_argument("--max-bin", type=int, default=127)
    ap.add_argument("--rounds", type=int, default=1000)
    ap.add_argument("--early-stop", type=int, default=50)
    ap.add_argument("--out", default=None)
    ap.add_argument("--plant-beta", type=float, default=0.0,
                    help="replace labels with ones generated from a known effect "
                         "of this size on a single real feature, to check the "
                         "model can find a signal that is actually there")
    ap.add_argument("--plant-feature", default="geo.moon.syn.h2.cos")
    a = ap.parse_args()
    out = a.out or (a.data + ".gbt.json")

    names, rows, X = cl.load(a.data)
    d = len(names)
    case_day = {s: dd for s, dd, c in zip(rows["stratum"], rows["day"], rows["case"]) if c}
    dos = np.array([case_day.get(s, np.nan) for s in rows["stratum"]])
    train = cl.Split(X, rows, dos < TRAIN_END, d, a.max_controls)
    val = cl.Split(X, rows, (dos >= TRAIN_END) & (dos < VAL_END), d, a.max_controls)
    print(f"train {train.n_strata} strata, validate {val.n_strata} strata")

    t0 = time.time()
    train.load(); val.load()
    print(f"loaded {(train.Xm.nbytes + val.Xm.nbytes)/1e9:.2f} GB "
          f"({time.time()-t0:.0f}s)")

    if a.plant_beta > 0:
        k = names.index(a.plant_feature)
        rng = np.random.default_rng(20260822)
        for sp in (train, val):
            col = sp.Xm[:, k].astype(np.float64)
            col = (col - col.mean()) / max(col.std(), 1e-12)
            y = np.zeros(len(col))
            for i in range(len(sp.starts)):
                lo, hi = sp.starts[i], sp.starts[i] + sp.sizes[i]
                w = np.exp(a.plant_beta * col[lo:hi])
                y[lo + rng.choice(sp.sizes[i], p=w / w.sum())] = 1.0
            sp.y = y
        print(f"PLANTED beta={a.plant_beta} on '{a.plant_feature}' -- "
              f"labels are synthetic, this is a capability check")

    t0 = time.time()
    dtrain = lgb.Dataset(train.Xm, label=train.y, params={"max_bin": a.max_bin},
                         free_raw_data=True)
    dtrain.construct()
    dval = lgb.Dataset(val.Xm, label=val.y, reference=dtrain, free_raw_data=True)
    dval.construct()
    # LightGBM has copied the data into its binned form; the float32 originals are
    # the memory bottleneck on this machine and are no longer needed.
    train.Xm = None
    val.Xm = None
    gc.collect()
    print(f"binned to {a.max_bin} bins ({time.time()-t0:.0f}s)")

    obj = make_objective(train.starts, train.sizes, train.seg, train.y)
    ev = make_eval([
        (train.starts, train.sizes, train.seg, train.y),
        (val.starts, val.sizes, val.seg, val.y),
    ])

    results = []
    # The pre-registered grid crosses trees x depth x learning rate. The tree count
    # is resolved by early stopping, which docs/22-model-grid.md permits, so the
    # sweep here is over depth and learning rate.
    for depth in (2, 3, 4):
        for lr in (0.01, 0.05):
            t0 = time.time()
            params = {
                "objective": obj,
                "learning_rate": lr,
                "max_depth": depth,
                "num_leaves": 2 ** depth,
                "min_data_in_leaf": 20,
                "bagging_fraction": 0.8,
                "bagging_freq": 1,
                "feature_fraction": 0.5,
                "max_bin": a.max_bin,
                "verbose": -1,
                "seed": 20260822,
                "num_threads": 6,
            }
            hist = {}
            booster = lgb.train(
                params, dtrain, num_boost_round=a.rounds,
                valid_sets=[dtrain, dval], valid_names=["train", "valid"],
                feval=ev,
                callbacks=[
                    lgb.early_stopping(a.early_stop, first_metric_only=False,
                                       verbose=False),
                    lgb.record_evaluation(hist),
                ],
            )
            best = booster.best_iteration or a.rounds
            ig_tr = hist["train"]["ig"][best - 1]
            ig_va = hist["valid"]["ig"][best - 1]
            results.append({"depth": depth, "lr": lr, "best_iter": int(best),
                            "train_ig": float(ig_tr), "val_ig": float(ig_va)})
            print(f"  depth {depth}  lr {lr:<5}  trees {best:>4}  "
                  f"train {ig_tr:+.5f}  validate {ig_va:+.5f} bits  "
                  f"({time.time()-t0:.0f}s)", flush=True)
            booster.save_model(a.data + f".gbt_d{depth}_lr{lr}.txt",
                               num_iteration=best)

    best = max(results, key=lambda r: r["val_ig"])
    med = float(np.median([r["val_ig"] for r in results]))
    print(f"\nbest on validation: depth {best['depth']}, lr {best['lr']}, "
          f"{best['val_ig']:+.5f} bits/event")
    print(f"median across {len(results)} configs: {med:+.5f} bits/event")
    json.dump({"data": a.data, "results": results, "best": best, "median": med},
              open(out, "w"), indent=2)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
