"""Stage 1: fixed predictors for the Rotterdam -> GBSG application (see SPECIFICATION.md).

Writes calibration inputs, target inputs containing only (q0, q1), and target outcomes in a separate
folder that is read only by stage 3. Patient-level files are written under work/rotterdam_gbsg/.
"""
import json
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
from common import (SEEDS, CAP, TAU, dump_json, feature_frame, fit_predictor,  # noqa: E402
                    predict_survival, serialize_lock, sha256, source_outcome, target_outcome)

DATA = ROOT / "data"
W = ROOT / "work" / "rotterdam_gbsg"
OUT = ROOT / "results" / "rotterdam_gbsg"
Q0 = ["age3", "age3log", "meno", "size_mid", "size_large", "hormon"]
Q1 = Q0 + ["nodes_inv_sqrt", "er"]
RAW_HASH = {
    "rotterdam.csv": "62703670d3be5d49c5476fb7131b020bdca9f1c25bb4522fd6abe8f9dab61f1b",
    "gbsg.csv": "9fa0fef0575d04d3273869be28405307f1b02993751498e8a0c844541cccd5d9",
}


def main():
    for name, expected in RAW_HASH.items():
        if sha256(DATA / name) != expected:
            raise RuntimeError("SHA-256 mismatch for data/" + name + "; see data/README.md")
    source = pd.read_csv(DATA / "rotterdam.csv")
    target = pd.read_csv(DATA / "gbsg.csv")
    source = source.loc[source["nodes"] > 0].reset_index(drop=True)
    stime, sevent = source_outcome(source)
    ttime, tevent = target_outcome(target)
    sf = feature_frame(source, "rotterdam")
    tf = feature_frame(target, "gbsg")
    analysis_root = W / "analysis_inputs"
    sealed_root = W / "sealed_target_outcomes"
    model_root = W / "model_locks"
    for directory in [analysis_root, sealed_root, model_root, OUT]:
        directory.mkdir(parents=True, exist_ok=True)
    sealed_path = sealed_root / "target_GBSG_outcomes.csv"
    pd.DataFrame({"patient_id": target["pid"].astype(str), "time": ttime, "event": tevent}).to_csv(sealed_path, index=False)
    outputs = []
    for seed in SEEDS:
        idx = np.arange(len(source))
        train_idx, cal_idx = train_test_split(idx, test_size=1/3, random_state=seed, stratify=sevent)
        lock0 = fit_predictor(sf.iloc[train_idx], stime[train_idx], sevent[train_idx], Q0, penalizer=0.01)
        lock1 = fit_predictor(sf.iloc[train_idx], stime[train_idx], sevent[train_idx], Q1, penalizer=0.01)
        seed_dir = analysis_root / ("seed{}".format(seed))
        seed_dir.mkdir(exist_ok=True)
        cal_path = seed_dir / "calibration_GBSG.csv"
        target_path = seed_dir / "target_GBSG.csv"
        pd.DataFrame({
            "patient_id": source.iloc[cal_idx]["pid"].astype(str).to_numpy(),
            "time": stime[cal_idx], "event": sevent[cal_idx],
            "q0_surv5y": predict_survival(lock0, sf.iloc[cal_idx]),
            "q1_surv5y": predict_survival(lock1, sf.iloc[cal_idx]),
        }).to_csv(cal_path, index=False)
        target_frame = pd.DataFrame({
            "patient_id": target["pid"].astype(str),
            "q0_surv5y": predict_survival(lock0, tf),
            "q1_surv5y": predict_survival(lock1, tf),
        })
        target_frame.to_csv(target_path, index=False)
        if set(target_frame.columns) != {"patient_id", "q0_surv5y", "q1_surv5y"}:
            raise RuntimeError("target inputs must not contain outcomes")
        model_path = model_root / ("seed{}.json".format(seed))
        dump_json({"seed": seed, "q0_columns": Q0, "q1_columns": Q1,
                   "n_train": int(len(train_idx)), "n_cal": int(len(cal_idx)),
                   "q0": serialize_lock(lock0), "q1": serialize_lock(lock1),
                   "tau_days": TAU, "cap_days": CAP}, model_path)
        outputs += [cal_path, target_path, model_path]
    manifest = {
        "status": "PASS",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "application": "Rotterdam -> GBSG",
        "target_outcomes_previously_examined_for_a_different_predictor_contrast": True,
        "target_outcomes_in_target_inputs": False,
        "specification_sha256": sha256(HERE / "SPECIFICATION.md"),
        "runner_sha256": sha256(__file__), "common_sha256": sha256(HERE / "common.py"),
        "raw_sha256": {k: sha256(DATA / k) for k in RAW_HASH},
        "analysis_outputs": [{"path": str(p.relative_to(ROOT)), "sha256": sha256(p)} for p in outputs],
        "sealed_target_outcomes": {"path": str(sealed_path.relative_to(ROOT)), "sha256": sha256(sealed_path)},
        "python": sys.version, "platform": platform.platform(),
    }
    dump_json(manifest, OUT / "STAGE1_MANIFEST.json")
    print(json.dumps({"status": manifest["status"], "source_n": len(source), "target_n": len(target)}, indent=2))


if __name__ == "__main__":
    main()
