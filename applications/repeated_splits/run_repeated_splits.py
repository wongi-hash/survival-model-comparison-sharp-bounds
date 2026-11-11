"""Repeated-split sensitivity analysis (see PROTOCOL.md). Resumable.

Usage:
  python run_repeated_splits.py run <app> [max_seconds]     app in {gbsg, temporal, nwts}
  python run_repeated_splits.py hash                        (record SHA-256 of the outcome-free records)
  python run_repeated_splits.py retro <app> [max_seconds]   (only after 'hash')
  python run_repeated_splits.py hash-retro
"""
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
import direct_vs_separate as D  # noqa: E402
from survival_tools import DATA, RESULTS, fit_ridge_cox, ipcw_brier_contrast, predict, sha256, stratified_split  # noqa: E402

OUT = RESULTS / "repeated_splits"
SPLITS = range(1001, 1101); BOOT = 999
LAB = {"ADOPT": "CAND", "KEEP_REFERENCE": "REF", "DEFER": "DEFER"}


def rot_features(raw, cohort):
    age10 = raw["age"].to_numpy(float) / 10.0
    if cohort == "rotterdam":
        size = raw["size"].astype(str); mid = size.eq("20-50").astype(float).to_numpy(); large = size.eq(">50").astype(float).to_numpy()
    else:
        s = raw["size"].to_numpy(float); mid = ((s > 20) & (s <= 50)).astype(float); large = (s > 50).astype(float)
    return pd.DataFrame({"age3": age10 ** 3, "age3log": age10 ** 3 * np.log(age10), "meno": raw["meno"].to_numpy(float),
                         "size_mid": mid, "size_large": large, "nodes_inv_sqrt": 1 / np.sqrt(raw["nodes"].to_numpy(float)),
                         "hormon": raw["hormon"].to_numpy(float), "er": raw["er"].to_numpy(float)})


def rot_rfs(raw, cap=365.25 * 7):
    t = np.where(raw["recur"].to_numpy() == 1, raw["rtime"].to_numpy(), raw["dtime"].to_numpy()).astype(float)
    e = ((raw["recur"].to_numpy() == 1) | (raw["death"].to_numpy() == 1)).astype(int)
    return np.minimum(t, cap), np.where(t <= cap, e, 0)


def load(app):
    RQ0 = ["age3", "age3log", "meno", "size_mid", "size_large", "hormon"]; RQ1 = RQ0 + ["nodes_inv_sqrt", "er"]
    if app in ("gbsg", "temporal"):
        rot = pd.read_csv(DATA / "rotterdam.csv"); rot = rot[rot.nodes > 0].reset_index(drop=True)
        if app == "gbsg":
            gb = pd.read_csv(DATA / "gbsg.csv")
            st, se = rot_rfs(rot); tt = np.minimum(gb.rfstime.to_numpy(float), 365.25 * 7)
            te = np.where(gb.rfstime.to_numpy(float) <= 365.25 * 7, gb.status.to_numpy(int), 0)
            return rot_features(rot, "rotterdam"), st, se, rot_features(gb, "gbsg"), tt, te, RQ0, RQ1, 365.25 * 5
        src = rot[rot.year <= 1988].reset_index(drop=True); tgt = rot[rot.year > 1988].reset_index(drop=True)
        st, se = rot_rfs(src); tt, te = rot_rfs(tgt)
        return rot_features(src, "rotterdam"), st, se, rot_features(tgt, "rotterdam"), tt, te, RQ0, RQ1, 365.25 * 5
    nw = pd.read_csv(DATA / "nwtco.csv"); cap = 365.25 * 5

    def f(r):
        s = r.stage.to_numpy(int)
        return pd.DataFrame({"stage2": (s == 2) * 1.0, "stage3": (s == 3) * 1.0, "stage4": (s == 4) * 1.0, "age_years": r.age.to_numpy(float) / 12,
                             "instit_unfav": (r.instit.to_numpy(int) == 2) * 1.0, "histol_unfav": (r.histol.to_numpy(int) == 2) * 1.0})

    def o(r):
        t = r.edrel.to_numpy(float); e = r.rel.to_numpy(int); return np.minimum(t, cap), np.where(t <= cap, e, 0)
    src = nw[nw.study == 3].reset_index(drop=True); tgt = nw[nw.study == 4].reset_index(drop=True)
    st, se = o(src); tt, te = o(tgt); Q0 = ["stage2", "stage3", "stage4", "age_years", "instit_unfav"]
    return f(src), st, se, f(tgt), tt, te, Q0, Q0 + ["histol_unfav"], 365.25 * 3


def predictors(app, s):
    SX, st, se, TX, tt, te, Q0, Q1, tau = load(app); dev, cal = stratified_split(se, s)
    q = {}
    for k, cols in (("q0", Q0), ("q1", Q1)):
        lock = fit_ridge_cox(SX.iloc[dev][cols].to_numpy(), st[dev], se[dev], tau)
        q[k] = (predict(lock, SX.iloc[cal][cols].to_numpy()), predict(lock, TX[cols].to_numpy()))
    calf = pd.DataFrame({"time": st[cal], "event": se[cal], "q0": q["q0"][0], "q1": q["q1"][0]})
    tgtf = pd.DataFrame({"q0": q["q0"][1], "q1": q["q1"][1]})
    return calf, tgtf, tt, te, tau


def run(app, limit=float("inf")):
    out = OUT / "records"; out.mkdir(parents=True, exist_ok=True); path = out / f"{app}.jsonl"
    done = {json.loads(l)["split"] for l in path.read_text().splitlines()} if path.exists() else set(); t0 = time.time()
    with path.open("a") as fh:
        for s in SPLITS:
            if s in done: continue
            if time.time() - t0 > limit: break
            rec = {"app": app, "split": s}
            try:
                calf, tgtf, _, _, tau = predictors(app, s)
                rec["n_cal"] = int(len(calf)); rec["cal_events_tau"] = int(((calf.time <= tau) & (calf.event == 1)).sum())
                rec["same_side"] = float((((tgtf.q0 - .5) * (tgtf.q1 - .5)) > 0).mean())
                rows, meta = D.analyze(app, "repsplit", s, calf, tgtf, tau, BOOT, rng_seed=D.stable_seed("repsplit-20261110", app, s, BOOT))
                rec["fit_failure"] = 0; rec["boot_valid"] = int(meta["bootstrap_valid"])
                rec["cells"] = [{k: (float(r[k]) if not isinstance(r[k], str) else LAB.get(r[k], r[k])) for k in
                                 ("gamma", "L_dir", "U_dir", "L_sep", "U_sep", "dir_sim_lower", "dir_sim_upper", "sep_sim_lower", "sep_sim_upper",
                                  "ID_width_ratio", "band_ratio", "action_dir", "action_sep")} for r in rows]
            except Exception as exc:
                rec["fit_failure"] = 1; rec["error"] = f"{type(exc).__name__}: {str(exc)[:200]}"
            fh.write(json.dumps(rec) + "\n"); fh.flush()
    print(app, sum(1 for _ in path.open()))


def retro(app, limit=float("inf")):
    if not (OUT / "OUTCOME_FREE_SHA256.txt").is_file():
        raise RuntimeError("run 'hash' first")
    path = OUT / "records" / f"{app}_retro.jsonl"
    done = {json.loads(l)["split"] for l in path.read_text().splitlines()} if path.exists() else set(); t0 = time.time()
    with path.open("a") as fh:
        for s in SPLITS:
            if s in done: continue
            if time.time() - t0 > limit: break
            calf, tgtf, tt, te, tau = predictors(app, s); q0, q1 = tgtf.q0.to_numpy(), tgtf.q1.to_numpy()
            est = ipcw_brier_contrast(tt, te, q0, q1, tau); rng = np.random.default_rng(20261110 + s); bs = []
            for _ in range(BOOT):
                i = rng.integers(0, len(tt), len(tt)); bs.append(ipcw_brier_contrast(tt[i], te[i], q0[i], q1[i], tau))
            lo, hi = np.quantile(bs, [0.025, 0.975])
            fh.write(json.dumps({"app": app, "split": s, "contrast": float(est), "ci_low": float(lo), "ci_high": float(hi)}) + "\n"); fh.flush()
    print(app, "retro", sum(1 for _ in path.open()))


def do_hash(retro_stage=False):
    files = sorted((OUT / "records").glob("*_retro.jsonl" if retro_stage else "*.jsonl"))
    if not retro_stage:
        files = [p for p in files if "retro" not in p.name]
    lines = [f"{sha256(p)}  {p.name}" for p in files] + [time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())]
    name = "RETRO_SHA256.txt" if retro_stage else "OUTCOME_FREE_SHA256.txt"
    (OUT / name).write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    cmd = sys.argv[1]
    lim = float(sys.argv[3]) if len(sys.argv) > 3 else float("inf")
    if cmd == "run": run(sys.argv[2], lim)
    elif cmd == "retro": retro(sys.argv[2], lim)
    elif cmd == "hash": do_hash()
    elif cmd == "hash-retro": do_hash(retro_stage=True)
    else: raise SystemExit(__doc__)
