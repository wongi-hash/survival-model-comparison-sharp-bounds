"""Simulation study: direct identification versus separate-risk subtraction (see PROTOCOL.md).

15 conditions: source calibration size n_S in {100, 200, 300, 400, 500} x overall source censoring in
{70%, 50%, 30%}; n_T = 500; 1,000 replications per condition; B = 199 bootstrap replicates.
Direct and separate endpoints use the same simulated data and bootstrap draws; each method has its own
maximum-statistic critical value. Records are appended to results/simulation/records/cell_<nS>_<pct>.jsonl
and the run is resumable.

Usage:
  python simulation/run_factorial.py <nS> <censoring_pct> [max_seconds]
  python simulation/run_factorial.py all
  python simulation/run_factorial.py summarize
"""
import json
import math
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import factorial_core as G  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "simulation"
GAM = G.GAMMAS; BOOT = G.BOOT; N_T = 500; REPS = 1000
NS = (100, 200, 300, 400, 500); PCT = (70, 50, 30)
expit, logit = G.expit, G.logit


def separate_endpoints(q0, q1, ps):
    """Separate-risk subtraction: bound each risk separately and subtract the intervals."""
    lp = logit(ps); out = np.empty((len(GAM), 2))
    for j, g in enumerate(GAM):
        pm, pp = expit(lp - g), expit(lp + g)

        def single(q):
            c = 1 - 2 * q
            return np.mean(q * q + c * np.where(c > 0, pm, pp)), np.mean(q * q + c * np.where(c > 0, pp, pm))
        l1, u1 = single(q1); l0, u0 = single(q0)
        out[j] = (l1 - u0, u1 - l0)
    return out


def truths():
    rng = np.random.default_rng(990011); x = rng.normal(G.TXSHIFT, 1, 2_000_000)
    q0, q1 = G.qv(x); ps = G.psrc(q0, q1)
    return G.endpoints(q0, q1, ps), separate_endpoints(q0, q1, ps)


T_DIR, T_SEP = truths()
assert np.allclose(T_DIR, G.TRUE_EP)


def decide(lo, up):
    return np.where(up < 0, "CAND", np.where(lo > 0, "REF", "DEFER"))


def one_rep(seed, nS, lam_c0):
    rng = np.random.default_rng(seed)
    xs = rng.normal(0, 1, nS); xt = rng.normal(G.TXSHIFT, 1, N_T)
    q0s, q1s = G.qv(xs); q0t, q1t = G.qv(xt); ps = G.psrc(q0s, q1s)
    lam = -np.log(np.clip(ps, 1e-8, 1)) / G.TAU ** G.WK
    fail = (-np.log(rng.random(nS)) / np.clip(lam, 1e-10, None)) ** (1 / G.WK)
    z = np.column_stack([q0s, q1s]); zs = (z - z.mean(0)) / np.where(z.std(0) < 1e-8, 1, z.std(0))
    cens = -np.log(rng.random(nS)) / (lam_c0 * np.exp(zs @ G.ALPHA_C))
    tm = np.minimum(fail, cens); ev = (fail <= cens).astype(int)
    rec = {"cens": float(1 - ev.mean()), "events": int(ev.sum()), "fatal": 0}
    try:
        psh = G.fitps(tm, ev, q0s, q1s, q0t, q1t)
    except Exception:
        rec["fatal"] = 1; return rec
    pt = {"dir": G.endpoints(q0t, q1t, psh), "sep": separate_endpoints(q0t, q1t, psh)}
    draws = {k: np.empty((BOOT, len(GAM), 2)) for k in pt}; v = 0
    rb = np.random.default_rng(seed + 800_000_003)
    for _ in range(BOOT):
        si = rb.integers(0, nS, nS); ti = rb.integers(0, N_T, N_T)
        try:
            psb = G.fitps(tm[si], ev[si], q0s[si], q1s[si], q0t, q1t)
            draws["dir"][v] = G.endpoints(q0t[ti], q1t[ti], psb[ti])
            draws["sep"][v] = separate_endpoints(q0t[ti], q1t[ti], psb[ti]); v += 1
        except Exception:
            pass
    if v < int(0.9 * BOOT):
        rec["fatal"] = 1; return rec
    for k, truth in (("dir", T_DIR), ("sep", T_SEP)):
        d = draws[k][:v]; se = np.maximum(d.std(0, ddof=1), 1e-10)
        crit = float(np.quantile(np.abs((d - pt[k][None]) / se[None]).max(axis=(1, 2)), 0.95))
        lo = pt[k][:, 0] - crit * se[:, 0]; up = pt[k][:, 1] + crit * se[:, 1]
        rec[f"{k}_width"] = (pt[k][:, 1] - pt[k][:, 0]).tolist()
        rec[f"{k}_outer"] = int(np.all((lo <= truth[:, 0]) & (truth[:, 1] <= up)))
        rec[f"{k}_decision"] = decide(lo, up).tolist(); rec[f"{k}_crit"] = crit
    return rec


def run(nS, pct, limit=float("inf")):
    out = OUT / "records"; out.mkdir(parents=True, exist_ok=True); path = out / f"cell_{nS}_{pct}.jsonl"
    done = {json.loads(l)["rep"] for l in path.read_text().splitlines()} if path.exists() else set()
    lam = G.calib_lambda(pct / 100); t0 = time.time()
    with path.open("a") as fh:
        for i in range(REPS):
            if i in done: continue
            if time.time() - t0 > limit: break
            r = one_rep(G.stable_seed("factorial-v2", nS, pct, i), nS, lam); r["rep"] = i
            fh.write(json.dumps(r) + "\n"); fh.flush()
    print(nS, pct, sum(1 for _ in path.open()))


def wilson(k, n, z=1.959963984540054):
    p = k / n; d = 1 + z * z / n; c = (p + z * z / (2 * n)) / d; h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return c - h, c + h


def summarize():
    import csv
    rows = []
    for nS in NS:
        for pct in PCT:
            R = [json.loads(l) for l in (OUT / "records" / f"cell_{nS}_{pct}.jsonl").read_text().splitlines()]
            ok = [r for r in R if not r["fatal"]]; n = len(ok)
            row = {"nS": nS, "cens_target": pct, "reps": len(R), "valid": n,
                   "mean_cens": round(np.mean([r["cens"] for r in ok]), 3), "mean_events": round(np.mean([r["events"] for r in ok]), 1)}
            for k in ("dir", "sep"):
                c = sum(r[f"{k}_outer"] for r in ok); lo, hi = wilson(c, n)
                row[f"{k}_coverage"] = round(c / n, 3); row[f"{k}_cov_lo"] = round(lo, 3); row[f"{k}_cov_hi"] = round(hi, 3)
            W = np.array([[r["sep_width"][j] / r["dir_width"][j] for j in range(1, len(GAM))] for r in ok])
            row["width_ratio_mean_g>0"] = round(float(W.mean()), 2)
            for j, gm in enumerate(GAM):
                for k in ("dir", "sep"):
                    dec = [r[f"{k}_decision"][j] for r in ok]
                    row[f"{k}_correct_g{gm}"] = round(dec.count("REF") / n, 3)
                    row[f"{k}_wrong_g{gm}"] = round(dec.count("CAND") / n, 4)
            rows.append(row)
    with open(OUT / "simulation_summary.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    print(f"[written] {OUT / 'simulation_summary.csv'}")


if __name__ == "__main__":
    a = sys.argv[1]
    if a == "summarize":
        summarize()
    elif a == "all":
        for nS in NS:
            for pct in PCT:
                run(nS, pct)
    else:
        run(int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3]) if len(sys.argv) > 3 else float("inf"))
