#!/usr/bin/env python3
"""Model B on a catalogue larger than memory, via streaming dataset construction.

    ~/eqf-work/venv/bin/python scripts/train_gbt_stream.py \
        --data /Volumes/2TB_EXT_1B/eqf-data/m40_timestrat --max-strata 45000

The M4.0 design matrix is 27 GB of float32. LightGBM's binned form is one byte
per value, which is 6.7 GB -- still beyond this machine. Two things make it fit:

  lgb.Sequence   builds the binned dataset by reading batches from the memmap,
                 so the float array is never materialised whole. Without it the
                 construction peak is the full 27 GB regardless of how small the
                 binned result is.
  --max-strata   subsamples whole strata. Rows within a stratum must stay
                 together: the objective is a softmax over the stratum, and a
                 stratum split across the sample would be scored against a
                 partial denominator.

Subsampling strata costs power in proportion, but 45,000 strata is still nearly
four times the entire M5.5+ catalogue.
"""

import argparse, gc, json, os, sys, time
import numpy as np
import lightgbm as lgb

TRAIN_END, VAL_END = 3652.0, 6210.0
BATCH = 4096


class MemmapSeq(lgb.Sequence):
    """Rows of one split, served to LightGBM in batches straight from the memmap."""

    def __init__(self, X, idx, batch_size=BATCH):
        self.X = X
        self.idx = idx
        self.batch_size = batch_size

    def __getitem__(self, k):
        if isinstance(k, slice):
            return np.asarray(self.X[self.idx[k]], dtype=np.float64)
        return np.asarray(self.X[self.idx[k]], dtype=np.float64)

    def __len__(self):
        return len(self.idx)


def strata_of(rows, n_all):
    stratum = rows["stratum"]
    starts = np.concatenate(([0], np.flatnonzero(np.diff(stratum)) + 1))
    sizes = np.diff(np.append(starts, n_all))
    return starts, sizes


def make_objective(starts, sizes, seg, y):
    def obj(preds, _d):
        m = np.maximum.reduceat(preds, starts)
        e = np.exp(preds - m[seg])
        p = e / np.add.reduceat(e, starts)[seg]
        return p - y, np.maximum(p * (1.0 - p), 1e-6)
    return obj


def make_eval(splits):
    by_len = {}
    for starts, sizes, seg, y in splits:
        n = int(np.sum(sizes))
        assert n not in by_len, "splits have equal row counts; cannot dispatch"
        by_len[n] = (starts, sizes, seg, y, -np.log(sizes.astype(np.float64)))

    def ev(preds, _d):
        starts, sizes, seg, y, null = by_len[len(preds)]
        m = np.maximum.reduceat(preds, starts)
        e = np.exp(preds - m[seg])
        lp = preds[y == 1] - m - np.log(np.add.reduceat(e, starts))
        return "ig", float(np.mean(lp - null) / np.log(2)), True
    return ev


def build_split(X, rows, sel_strata, starts, sizes, case, d):
    """Row indices, stratum layout and labels for a chosen set of strata."""
    idx = np.concatenate([np.arange(starts[s], starts[s] + sizes[s]) for s in sel_strata])
    sz = sizes[sel_strata]
    st = np.concatenate(([0], np.cumsum(sz)[:-1]))
    seg = np.repeat(np.arange(len(sz)), sz)
    y = case[idx].astype(np.float64)
    return idx, st, sz, seg, y


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True)
    ap.add_argument("--max-strata", type=int, default=45000)
    ap.add_argument("--max-bin", type=int, default=63)
    ap.add_argument("--rounds", type=int, default=1000)
    ap.add_argument("--early-stop", type=int, default=50)
    ap.add_argument("--plant-beta", type=float, default=0.0)
    ap.add_argument("--plant-sweep", default=None,
                    help="comma-separated betas; bins the data once and refits "
                         "for each, since only the labels change")
    ap.add_argument("--plant-feature", default="geo.moon.syn.h2.cos")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    out = a.out or (a.data + ".gbt.json")

    names = [l.strip() for l in open(a.data + ".names") if l.strip()]
    d = len(names)
    rows = np.genfromtxt(a.data + ".rows.csv", delimiter=",", names=True,
                         dtype=None, encoding="utf-8")
    n_all = len(rows)
    X = np.memmap(a.data + ".f32", dtype=np.float32, mode="r", shape=(n_all, d))
    starts, sizes = strata_of(rows, n_all)
    case = rows["case"].astype(bool)
    day = rows["day"]
    print(f"{n_all} rows, {len(starts)} strata, {d} features")

    case_day = np.full(len(starts), np.nan)
    for si, (lo, sz) in enumerate(zip(starts, sizes)):
        c = np.flatnonzero(case[lo:lo + sz])
        if len(c) == 1:
            case_day[si] = day[lo + c[0]]
    usable = (sizes > 1) & ~np.isnan(case_day)
    tr_s = np.flatnonzero(usable & (case_day < TRAIN_END))
    va_s = np.flatnonzero(usable & (case_day >= TRAIN_END) & (case_day < VAL_END))

    rng = np.random.default_rng(20260822)
    if len(tr_s) > a.max_strata:
        keep = rng.choice(len(tr_s), a.max_strata, replace=False)
        keep.sort()
        tr_s = tr_s[keep]
        print(f"subsampled training strata to {len(tr_s)}")
    print(f"train {len(tr_s)} strata, validate {len(va_s)} strata (test sealed)")

    tr_idx, tr_st, tr_sz, tr_seg, tr_y = build_split(X, rows, tr_s, starts, sizes, case, d)
    va_idx, va_st, va_sz, va_seg, va_y = build_split(X, rows, va_s, starts, sizes, case, d)

    def replant(beta):
        """Redraw both splits' labels from a known effect on one real feature.

        The feature matrix is untouched, so the binned Dataset built below stays
        valid across the whole sweep -- with a custom objective LightGBM reads the
        label only through the closures, never from the Dataset itself.
        """
        k = names.index(a.plant_feature)
        for idx, st, sz, y in ((tr_idx, tr_st, tr_sz, tr_y), (va_idx, va_st, va_sz, va_y)):
            col = np.asarray(X[idx, k], dtype=np.float64)
            col = (col - col.mean()) / max(col.std(), 1e-12)
            y[:] = 0.0
            for i in range(len(st)):
                lo, hi = st[i], st[i] + sz[i]
                w = np.exp(beta * col[lo:hi])
                y[lo + rng.choice(sz[i], p=w / w.sum())] = 1.0

    if a.plant_beta > 0:
        print(f"PLANTED beta={a.plant_beta} on '{a.plant_feature}' -- "
              f"synthetic labels, this is a capability check")
        replant(a.plant_beta)

    t0 = time.time()
    params = {"max_bin": a.max_bin, "min_data_in_bin": 1, "verbose": -1}
    dtrain = lgb.Dataset(MemmapSeq(X, tr_idx), label=tr_y, params=params,
                         free_raw_data=True)
    dtrain.construct()
    dval = lgb.Dataset(MemmapSeq(X, va_idx), label=va_y, reference=dtrain,
                       params=params, free_raw_data=True)
    dval.construct()
    gc.collect()
    print(f"binned to {a.max_bin} bins ({time.time()-t0:.0f}s)")

    obj = make_objective(tr_st, tr_sz, tr_seg, tr_y)
    ev = make_eval([(tr_st, tr_sz, tr_seg, tr_y), (va_st, va_sz, va_seg, va_y)])

    def fit(depth, lr, rounds=None):
        p = {"objective": obj, "learning_rate": lr, "max_depth": depth,
             "num_leaves": 2 ** depth, "min_data_in_leaf": 20,
             "bagging_fraction": 0.8, "bagging_freq": 1,
             "feature_fraction": 0.5, "max_bin": a.max_bin,
             "verbose": -1, "seed": 20260822, "num_threads": 6}
        hist = {}
        bst = lgb.train(p, dtrain, num_boost_round=rounds or a.rounds,
                        valid_sets=[dtrain, dval], valid_names=["train", "valid"],
                        feval=ev,
                        callbacks=[lgb.early_stopping(a.early_stop, verbose=False),
                                   lgb.record_evaluation(hist)])
        b = bst.best_iteration or (rounds or a.rounds)
        return b, hist["train"]["ig"][b - 1], hist["valid"]["ig"][b - 1]

    if a.plant_sweep:
        # The power curve for Model B. Its null on real data means nothing without
        # knowing what size of effect it would have caught.
        print(f"\n{'beta':>6}  {'trees':>6} {'train IG':>10} {'valid IG':>10}   verdict")
        sweep = []
        for b in [float(x) for x in a.plant_sweep.split(",")]:
            replant(b) if b > 0 else replant(0.0)
            n, ig_tr, ig_va = fit(3, 0.05)
            verdict = ("DETECTED" if ig_va >= 0.01 else
                       "weak" if ig_va > 0 else "invisible")
            print(f"{b:>6.2f}  {n:>6} {ig_tr:>+10.5f} {ig_va:>+10.5f}   {verdict}",
                  flush=True)
            sweep.append({"beta": b, "trees": int(n), "train_ig": float(ig_tr),
                          "val_ig": float(ig_va)})
        json.dump({"data": a.data, "sweep": sweep,
                   "n_train_strata": int(len(tr_s))}, open(out, "w"), indent=2)
        print(f"wrote {out}")
        return

    results = []
    for depth in (2, 3, 4):
        for lr in (0.01, 0.05):
            t0 = time.time()
            best, ig_tr, ig_va = fit(depth, lr)
            results.append({"depth": depth, "lr": lr, "best_iter": int(best),
                            "train_ig": float(ig_tr), "val_ig": float(ig_va)})
            print(f"  depth {depth}  lr {lr:<5}  trees {best:>4}  "
                  f"train {ig_tr:+.5f}  validate {ig_va:+.5f} bits  "
                  f"({time.time()-t0:.0f}s)", flush=True)

    best = max(results, key=lambda r: r["val_ig"])
    med = float(np.median([r["val_ig"] for r in results]))
    print(f"\nbest on validation: depth {best['depth']}, lr {best['lr']}, "
          f"{best['val_ig']:+.5f} bits/event")
    print(f"median across {len(results)} configs: {med:+.5f} bits/event")
    json.dump({"data": a.data, "results": results, "best": best, "median": med,
               "n_train_strata": int(len(tr_s)), "n_val_strata": int(len(va_s)),
               "planted": a.plant_beta},
              open(out, "w"), indent=2)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
