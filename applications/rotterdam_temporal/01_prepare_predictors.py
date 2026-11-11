"""Stage 1: fixed predictors for the Rotterdam earlier -> later application (see PROTOCOL.md).

Source: node-positive Rotterdam patients with surgery in 1978-1988; target: surgery in 1989-1993.
Writes calibration inputs, target inputs containing only (q0, q1), and target outcomes in a separate
folder that is read only by stage 3. Patient-level files are written under work/rotterdam_temporal/.
"""
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from survival_tools import DATA, RESULTS, WORK, check_sha256, fit_ridge_cox, predict, sha256, stratified_split  # noqa: E402

HERE = Path(__file__).resolve().parent
RAW = DATA / "rotterdam.csv"
RAW_SHA = "62703670d3be5d49c5476fb7131b020bdca9f1c25bb4522fd6abe8f9dab61f1b"
TAU = 365.25 * 5
CAP = 365.25 * 7
SEEDS = (42, 52, 62, 72)
SPLIT_YEAR = 1988
Q0 = ["age3", "age3log", "meno", "size_mid", "size_large", "hormon"]
Q1 = Q0 + ["nodes_inv_sqrt", "er"]
W = WORK / "rotterdam_temporal"
OUT = RESULTS / "rotterdam_temporal"


def rfs(raw):
    """Recurrence-free survival (recurrence or death), administratively censored at seven years."""
    time = np.where(raw["recur"].to_numpy() == 1, raw["rtime"].to_numpy(), raw["dtime"].to_numpy()).astype(float)
    event = ((raw["recur"].to_numpy() == 1) | (raw["death"].to_numpy() == 1)).astype(int)
    event = np.where(time <= CAP, event, 0)
    return np.minimum(time, CAP), event


def features(raw):
    age10 = raw["age"].to_numpy(float) / 10.0
    size = raw["size"].astype(str)
    nodes = raw["nodes"].to_numpy(float)
    assert np.all(nodes > 0)
    return pd.DataFrame({
        "age3": age10 ** 3, "age3log": (age10 ** 3) * np.log(age10), "meno": raw["meno"].to_numpy(float),
        "size_mid": size.eq("20-50").astype(float).to_numpy(), "size_large": size.eq(">50").astype(float).to_numpy(),
        "nodes_inv_sqrt": 1.0 / np.sqrt(nodes), "hormon": raw["hormon"].to_numpy(float), "er": raw["er"].to_numpy(float),
    })


def main():
    check_sha256(RAW, RAW_SHA)
    raw = pd.read_csv(RAW)
    raw = raw.loc[raw["nodes"] > 0].reset_index(drop=True)
    src = raw.loc[raw["year"] <= SPLIT_YEAR].reset_index(drop=True)
    tgt = raw.loc[raw["year"] > SPLIT_YEAR].reset_index(drop=True)
    stime, sevent = rfs(src)
    sf, tf = features(src), features(tgt)
    inputs = W / "analysis_inputs"; sealed = W / "sealed_target_outcomes"
    for d in (inputs, sealed, OUT):
        d.mkdir(parents=True, exist_ok=True)
    ttime, tevent = rfs(tgt)
    sealed_path = sealed / "target_outcomes.csv"
    pd.DataFrame({"pid": tgt["pid"], "time": ttime, "event": tevent}).to_csv(sealed_path, index=False)
    manifest = {"created_utc": datetime.now(timezone.utc).isoformat(), "raw_sha256": RAW_SHA,
                "source_n": int(len(src)), "target_n": int(len(tgt)), "source_events": int(sevent.sum()),
                "protocol_sha256": sha256(HERE / "PROTOCOL.md"), "runner_sha256": sha256(__file__),
                "sealed_target_outcomes_sha256": sha256(sealed_path), "seeds": {}}
    for seed in SEEDS:
        dev, cal = stratified_split(sevent, seed)
        l0 = fit_ridge_cox(sf.iloc[dev][Q0].to_numpy(), stime[dev], sevent[dev], TAU)
        l1 = fit_ridge_cox(sf.iloc[dev][Q1].to_numpy(), stime[dev], sevent[dev], TAU)
        d = inputs / f"seed{seed}"; d.mkdir(exist_ok=True)
        calf = pd.DataFrame({"time": stime[cal], "event": sevent[cal],
                             "q0": predict(l0, sf.iloc[cal][Q0].to_numpy()), "q1": predict(l1, sf.iloc[cal][Q1].to_numpy())})
        tarf = pd.DataFrame({"q0": predict(l0, tf[Q0].to_numpy()), "q1": predict(l1, tf[Q1].to_numpy())})
        assert list(tarf.columns) == ["q0", "q1"]
        calf.to_csv(d / "calibration.csv", index=False); tarf.to_csv(d / "target.csv", index=False)
        manifest["seeds"][str(seed)] = {"n_dev": int(len(dev)), "n_cal": int(len(cal)), "cal_events": int(sevent[cal].sum()),
                                        "q0_beta": l0["beta"].round(6).tolist(), "q1_beta": l1["beta"].round(6).tolist(),
                                        "calibration_sha256": sha256(d / "calibration.csv"), "target_sha256": sha256(d / "target.csv")}
    (OUT / "STAGE1_MANIFEST.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(json.dumps({k: manifest[k] for k in ("source_n", "target_n", "source_events")}),
          {s: (v["n_cal"], v["cal_events"]) for s, v in manifest["seeds"].items()})


if __name__ == "__main__":
    main()
