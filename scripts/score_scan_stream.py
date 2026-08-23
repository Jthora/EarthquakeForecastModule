#!/usr/bin/env python3
"""Streaming conditional score test — full power, on catalogues too large to hold.

    ~/eqf-work/venv/bin/python scripts/score_scan_stream.py \
        --data /Volumes/2TB_EXT_1B/eqf-data/m40_dayoffset --perms 200

The M4.0 design matrix is 26.7 GB against 2.5 GB of usable RAM, and it is the
matrix that matters: 154,302 independent events instead of 12,160. The score
test does not need the data resident, because every quantity it needs is local
to a stratum:

    U_k      = sum_s ( x[case] - mean_s(x) )
    V_k      = sum_s mean_s( (x - mean_s(x))^2 )
    U_perm_k = sum_s ( x[random row of s] - mean_s(x) )

so a single pass over chunks aligned to stratum boundaries accumulates all of
them. Memory is one chunk plus the accumulators, regardless of catalogue size.

The permutations are the part that would normally force many passes. They do not
here: with W the P x n indicator matrix picking one row per stratum per
permutation, every permutation's statistic is a single sparse-dense product
U_perm = W @ Xc, so all P of them ride along in the same pass. 200 permutations
cost about 15 s of arithmetic on top of an I/O-bound 6-minute read.

Permutation rather than Bonferroni because cos and sin of neighbouring harmonics
of the same pair ask nearly the same question; permuting which row of a stratum
is the case preserves that correlation exactly.

The test period is never read.
"""

import argparse, json, os, sys, time
import numpy as np

TEST_START = 6210.0     # 2017-01-01
TARGET_CHUNK = 8192


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True)
    ap.add_argument("--perms", type=int, default=200)
    ap.add_argument("--top", type=int, default=30)
    ap.add_argument("--out", default=None)
    ap.add_argument("--plant-beta", type=float, default=0.0,
                    help="reassign which row of each stratum is the case, from a "
                         "known effect of this size on --plant-feature, to measure "
                         "what size of effect this scan would actually detect")
    ap.add_argument("--plant-feature", default="geo.moon.syn.h2.cos")
    a = ap.parse_args()
    out = a.out or (a.data + ".scan.json")

    names = [l.strip() for l in open(a.data + ".names") if l.strip()]
    d = len(names)
    rows = np.genfromtxt(a.data + ".rows.csv", delimiter=",", names=True,
                         dtype=None, encoding="utf-8")
    n_all = len(rows)
    X = np.memmap(a.data + ".f32", dtype=np.float32, mode="r", shape=(n_all, d))
    print(f"{n_all} rows x {d} features "
          f"({os.path.getsize(a.data + '.f32')/1e9:.1f} GB on disk)")

    stratum = rows["stratum"]
    case = rows["case"].astype(bool)
    day = rows["day"]

    # Stratum boundaries over the whole file, then keep only pre-test strata.
    starts = np.concatenate(([0], np.flatnonzero(np.diff(stratum)) + 1))
    sizes = np.diff(np.append(starts, n_all))
    case_day = np.full(len(starts), np.nan)
    for si, (lo, sz) in enumerate(zip(starts, sizes)):
        c = np.flatnonzero(case[lo:lo + sz])
        if len(c) == 1:
            case_day[si] = day[lo + c[0]]
    keep = (case_day < TEST_START) & (sizes > 1) & ~np.isnan(case_day)
    ks = starts[keep]
    kz = sizes[keep]
    n_strata = len(ks)
    print(f"{n_strata} usable strata before {TEST_START:.0f} "
          f"({int(kz.sum())} rows); test period never read")

    # Planting has to happen before the pass, and needs the chosen feature's column
    # standardised across the whole period, so it is read once up front.
    planted_case = None
    if a.plant_beta > 0:
        kf = names.index(a.plant_feature)
        col = np.empty(int(kz.sum()) + 0, dtype=np.float64)
        vals = []
        for lo, sz in zip(ks, kz):
            vals.append(np.asarray(X[lo:lo + sz, kf], dtype=np.float64))
        col = np.concatenate(vals)
        col = (col - col.mean()) / max(col.std(), 1e-12)
        prng = np.random.default_rng(20260823)
        planted_case = {}
        off = 0
        for si, (lo, sz) in enumerate(zip(ks, kz)):
            w = np.exp(a.plant_beta * col[off:off + sz])
            planted_case[si] = int(prng.choice(sz, p=w / w.sum()))
            off += sz
        print(f"PLANTED beta={a.plant_beta} on '{a.plant_feature}' -- "
              f"case rows reassigned; this measures detectable effect size")

    rng = np.random.default_rng(20260822)
    P = a.perms
    U = np.zeros(d, dtype=np.float64)
    V = np.zeros(d, dtype=np.float64)
    Up = np.zeros((P, d), dtype=np.float64)
    # Overall (not within-stratum) moments, to judge which features the design
    # actually has leverage on.
    T1 = np.zeros(d, dtype=np.float64)
    T2 = np.zeros(d, dtype=np.float64)
    n_rows_total = 0

    # Walk the strata in blocks whose rows are contiguous in the file, so each
    # read is sequential -- on a 76 MB/s external drive a scattered read would
    # dominate everything else.
    t0 = time.time()
    i = 0
    done_rows = 0
    while i < n_strata:
        j = i
        n_rows = 0
        while j < n_strata and n_rows < TARGET_CHUNK:
            # Stop the block if the next stratum is not adjacent in the file.
            if j > i and ks[j] != ks[j - 1] + kz[j - 1]:
                break
            n_rows += kz[j]
            j += 1
        lo, hi = ks[i], ks[j - 1] + kz[j - 1]
        # float64, not float32. Centring a slow feature in float32 is catastrophic:
        # an outer-planet aspect changes by less than float32 epsilon across a
        # 5-day stratum, so the centred residual is pure quantisation noise and
        # z = U/sqrt(V) explodes. Reading as f64 stops the arithmetic adding to
        # the problem; the filter below removes what the stored f32 already lost.
        blk = np.array(X[lo:hi], dtype=np.float64)

        loc_starts = (ks[i:j] - lo).astype(np.int64)
        loc_sizes = kz[i:j]
        seg = np.repeat(np.arange(j - i), loc_sizes)

        T1 += blk.sum(0)
        T2 += (blk * blk).sum(0)
        n_rows_total += hi - lo

        s = np.add.reduceat(blk, loc_starts, axis=0)
        means = s / loc_sizes[:, None]
        blk -= means[seg]                                  # centre in place

        sq = np.add.reduceat(blk * blk, loc_starts, axis=0)
        V += (sq / loc_sizes[:, None]).sum(0)

        if planted_case is None:
            c_rows = np.array([lst + np.flatnonzero(case[ks[i + t]:ks[i + t] + kz[i + t]])[0]
                               for t, lst in enumerate(loc_starts)])
        else:
            c_rows = np.array([lst + planted_case[i + t]
                               for t, lst in enumerate(loc_starts)])
        U += blk[c_rows].sum(0, dtype=np.float64)

        # One indicator row per permutation, then a single BLAS product carries
        # every permutation through this chunk at once.
        picks = loc_starts[None, :] + (rng.random((P, j - i)) * loc_sizes).astype(np.int64)
        W = np.zeros((P, hi - lo), dtype=np.float64)
        np.put_along_axis(W, picks, 1.0, axis=1)
        Up += W @ blk

        done_rows += hi - lo
        i = j
        if (i // 2000) != ((i - (j - i)) // 2000):
            el = time.time() - t0
            print(f"  {done_rows}/{int(kz.sum())} rows  {el:.0f}s elapsed, "
                  f"{el * (kz.sum()/max(done_rows,1) - 1):.0f}s left", flush=True)

    # A feature is testable only where the design moves it. Outer-planet aspects
    # barely change across a 5-day matched set; their within-stratum variance is
    # a rounding artefact of the stored f32, not information. Including them adds
    # nothing but inflates the multiplicity and the permutation threshold, which
    # is how a median permutation max |z| of 21 arises. The pre-registration
    # already states this design is blind to anything slower than its offsets --
    # this makes that blindness explicit instead of letting it pollute the test.
    total_var = np.maximum(T2 / n_rows_total - (T1 / n_rows_total) ** 2, 0.0)
    within_var = V / n_strata
    with np.errstate(divide="ignore", invalid="ignore"):
        leverage = np.where(total_var > 0, within_var / total_var, 0.0)
    ok = (V > 1e-12) & (leverage > 1e-8)
    n_dropped = int((V > 1e-12).sum() - ok.sum())
    print(f"{n_dropped} features dropped for no within-stratum leverage "
          f"(within-stratum sd below 1e-4 of overall sd)")
    z = np.zeros(d)
    z[ok] = U[ok] / np.sqrt(V[ok])
    zp = np.zeros((P, d))
    zp[:, ok] = Up[:, ok] / np.sqrt(V[ok])
    print(f"one pass in {time.time()-t0:.0f}s; {int(ok.sum())} testable features")

    zo = z[ok]
    print(f"\nz distribution: mean {zo.mean():+.4f}, sd {zo.std():.4f} "
          f"(expected 0.000, 1.000)")
    from scipy.stats import norm
    for t in (2, 3, 4):
        print(f"  |z| > {t}: {int((np.abs(zo) > t).sum())} of {int(ok.sum())} "
              f"({100*np.mean(np.abs(zo) > t):.3f}%, expect "
              f"{100*2*norm.sf(t):.3f}%)")

    obs = float(np.abs(zo).max())
    maxes = np.abs(zp[:, ok]).max(1)
    thresh = float(np.quantile(maxes, 0.95))
    pval = float((np.sum(maxes >= obs) + 1) / (P + 1))
    print(f"\npermutation null for max |z| over {P} draws:")
    print(f"  median {np.median(maxes):.3f}, 95th percentile {thresh:.3f}")
    print(f"  observed {obs:.3f}  ->  family-wise p = {pval:.4f}")
    # Permutation sd of z, as a check that the analytic variance is right.
    # The across-feature sd of a SINGLE draw, which is what the observed z is.
    # Features are heavily correlated, so one draw's spread is itself variable;
    # comparing the observed spread to the permutation draws' spreads says
    # whether an under-dispersed observed z field is unusual or ordinary.
    per_draw_sd = zp[:, ok].std(axis=1)
    print(f"  permutation z sd: {zp[:, ok].std():.4f} pooled; "
          f"per draw {per_draw_sd.mean():.3f} +/- {per_draw_sd.std():.3f}")
    sd_p = float((np.sum(per_draw_sd <= zo.std()) + 1) / (P + 1))
    print(f"  observed z sd {zo.std():.3f} sits at permutation quantile {sd_p:.3f}")

    order = np.argsort(-np.abs(z))
    print(f"\ntop {a.top} features by |z|:")
    for k in order[:a.top]:
        mark = " *" if abs(z[k]) >= thresh else ""
        print(f"  {z[k]:+8.3f}  {names[k]}{mark}")

    json.dump({"data": a.data, "n_strata": int(n_strata),
               "n_features": int(ok.sum()), "perms": P,
               "z_mean": float(zo.mean()), "z_sd": float(zo.std()),
               "perm_z_sd": float(zp[:, ok].std()),
               "max_abs_z": obs, "perm_threshold_95": thresh,
               "family_wise_p": pval,
               "top": [{"name": names[k], "z": float(z[k])} for k in order[:200]]},
              open(out, "w"), indent=2)
    print(f"\nwrote {out}")
    print("At least one feature separates cases from controls beyond chance."
          if pval < 0.05 else
          "No feature separates cases from controls beyond chance.")


if __name__ == "__main__":
    main()
