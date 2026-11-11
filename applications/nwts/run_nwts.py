"""NWTS-3 -> NWTS-4 application (see PROTOCOL.md).

Usage:
  python run_nwts.py prepare          # stage 1: fixed predictors; target outcomes stored separately
  python run_nwts.py evaluate <seed>  # stage 2: outcome-free direct vs separate (seed 42, 52, 62 or 72)
  python run_nwts.py aggregate        # stage 2: combine seeds and record SHA-256 hashes
  python run_nwts.py retrospective    # stage 3: only after stage 2 outputs are hashed
"""
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
import direct_vs_separate as D  # noqa: E402
from survival_tools import (DATA, RESULTS, WORK, check_sha256, fit_ridge_cox, ipcw_brier_scores,  # noqa: E402
                            predict, retrospective_interval, sha256, stratified_split)

HERE = Path(__file__).resolve().parent
RAW = DATA / "nwtco.csv"
RAW_SHA = "b5fa39f11d2adb8f51845dd302806da7ed09c4eecb453548a70f0692ebf77549"
TAU = 365.25 * 3
CAP = 365.25 * 5
SEEDS = (42, 52, 62, 72)
BOOT = 999
Q0 = ["stage2", "stage3", "stage4", "age_years", "instit_unfav"]
Q1 = Q0 + ["histol_unfav"]
LABELS = {"ADOPT": "CANDIDATE SUPPORTED", "KEEP_REFERENCE": "REFERENCE SUPPORTED", "DEFER": "DEFER"}
W = WORK / "nwts"
OUT = RESULTS / "nwts"


def outcome(raw):
    """Relapse-free survival, administratively censored at five years."""
    t = raw["edrel"].to_numpy(float); e = raw["rel"].to_numpy(int)
    e = np.where(t <= CAP, e, 0)
    return np.minimum(t, CAP), e


def features(raw):
    st = raw["stage"].to_numpy(int)
    return pd.DataFrame({"stage2": (st == 2).astype(float), "stage3": (st == 3).astype(float), "stage4": (st == 4).astype(float),
                         "age_years": raw["age"].to_numpy(float) / 12.0,
                         "instit_unfav": (raw["instit"].to_numpy(int) == 2).astype(float),
                         "histol_unfav": (raw["histol"].to_numpy(int) == 2).astype(float)})


def prepare():
    check_sha256(RAW, RAW_SHA)
    raw = pd.read_csv(RAW)
    src = raw[raw.study == 3].reset_index(drop=True); tgt = raw[raw.study == 4].reset_index(drop=True)
    stime, sevent = outcome(src); sf, tf = features(src), features(tgt)
    inp = W / "analysis_inputs"; sealed = W / "sealed_target_outcomes"
    for d in (inp, sealed, OUT):
        d.mkdir(parents=True, exist_ok=True)
    tt, te = outcome(tgt)
    pd.DataFrame({"seqno": tgt.seqno, "time": tt, "event": te}).to_csv(sealed / "target_outcomes.csv", index=False)
    man = {"created_utc": datetime.now(timezone.utc).isoformat(), "raw_sha256": RAW_SHA, "source_n": int(len(src)),
           "target_n": int(len(tgt)), "source_events": int(sevent.sum()), "protocol_sha256": sha256(HERE / "PROTOCOL.md"),
           "sealed_sha256": sha256(sealed / "target_outcomes.csv"), "seeds": {}}
    for seed in SEEDS:
        dev, cal = stratified_split(sevent, seed)
        l0 = fit_ridge_cox(sf.iloc[dev][Q0].to_numpy(), stime[dev], sevent[dev], TAU)
        l1 = fit_ridge_cox(sf.iloc[dev][Q1].to_numpy(), stime[dev], sevent[dev], TAU)
        d = inp / f"seed{seed}"; d.mkdir(exist_ok=True)
        pd.DataFrame({"time": stime[cal], "event": sevent[cal], "q0": predict(l0, sf.iloc[cal][Q0].to_numpy()),
                      "q1": predict(l1, sf.iloc[cal][Q1].to_numpy())}).to_csv(d / "calibration.csv", index=False)
        tarf = pd.DataFrame({"q0": predict(l0, tf[Q0].to_numpy()), "q1": predict(l1, tf[Q1].to_numpy())})
        assert list(tarf.columns) == ["q0", "q1"]; tarf.to_csv(d / "target.csv", index=False)
        man["seeds"][str(seed)] = {"n_dev": int(len(dev)), "n_cal": int(len(cal)), "cal_events": int(sevent[cal].sum()),
                                   "cal_events_before_tau": int(((stime[cal] <= TAU) & (sevent[cal] == 1)).sum()),
                                   "q0_beta": l0["beta"].round(6).tolist(), "q1_beta": l1["beta"].round(6).tolist()}
    (OUT / "STAGE1_MANIFEST.json").write_text(json.dumps(man, indent=2), encoding="utf-8")
    print(json.dumps({k: man[k] for k in ("source_n", "target_n", "source_events")}),
          {s: v["cal_events_before_tau"] for s, v in man["seeds"].items()})


def evaluate(seed):
    d = W / "analysis_inputs" / f"seed{seed}"
    cal = pd.read_csv(d / "calibration.csv"); target = pd.read_csv(d / "target.csv")
    if set(target.columns) != {"q0", "q1"}:
        raise RuntimeError("target inputs must contain only q0 and q1")
    rows, meta = D.analyze("NWTS3-NWTS4", "central_histology", seed, cal, target, TAU, BOOT,
                           rng_seed=D.stable_seed("nwtco-20261110", seed, BOOT))
    OUT.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(rows); df["decision_direct"] = df.action_dir.map(LABELS); df["decision_separate"] = df.action_sep.map(LABELS)
    df.to_csv(OUT / f"seed{seed}.csv", index=False)
    (OUT / f"seed{seed}_meta.json").write_text(json.dumps(meta, indent=2, default=str), encoding="utf-8")
    print(df[["seed", "gamma", "L_dir", "U_dir", "dir_sim_lower", "dir_sim_upper", "sep_sim_lower", "sep_sim_upper",
              "ID_width_ratio", "band_ratio", "decision_direct", "decision_separate"]].round(4).to_string(index=False))


def aggregate():
    df = pd.concat([pd.read_csv(OUT / f"seed{s}.csv") for s in SEEDS], ignore_index=True)
    df.to_csv(OUT / "outcome_free_all_cells.csv", index=False)
    files = [OUT / "outcome_free_all_cells.csv"] + [OUT / f"seed{s}.csv" for s in SEEDS]
    lines = [f"{sha256(p)}  {p.name}" for p in files] + [datetime.now(timezone.utc).isoformat()]
    (OUT / "OUTCOME_FREE_SHA256.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(df[["seed", "gamma", "decision_direct", "decision_separate"]].to_string(index=False))


def retrospective():
    if not (OUT / "OUTCOME_FREE_SHA256.txt").is_file():
        raise RuntimeError("run 'aggregate' first")
    oc = pd.read_csv(W / "sealed_target_outcomes" / "target_outcomes.csv")
    t, e = oc.time.to_numpy(float), oc.event.to_numpy(int)
    rows = []
    for seed in SEEDS:
        tg = pd.read_csv(W / "analysis_inputs" / f"seed{seed}" / "target.csv")
        q0, q1 = tg.q0.to_numpy(), tg.q1.to_numpy()
        b0, b1 = ipcw_brier_scores(t, e, q0, q1, TAU)
        est, lo, hi = retrospective_interval(t, e, q0, q1, TAU, seed=20261110 + seed, boot=BOOT)
        rows.append({"seed": seed, "brier_ref": b0, "brier_cand": b1, "contrast": est, "ci_low": lo, "ci_high": hi,
                     "target_censored_before_tau": float(((t <= TAU) & (e == 0)).mean())})
    df = pd.DataFrame(rows); df.to_csv(OUT / "retrospective_target_contrast.csv", index=False)
    print(df.round(4).to_string(index=False))


if __name__ == "__main__":
    a = sys.argv[1]
    {"prepare": prepare, "aggregate": aggregate, "retrospective": retrospective}.get(a, lambda: evaluate(int(sys.argv[2])))()
