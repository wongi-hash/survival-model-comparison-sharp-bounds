"""Stage 2: outcome-free direct versus separate evaluation (target inputs contain q0 and q1 only).

Usage: python 02_outcome_free_evaluation.py <seed>     (42, 52, 62 or 72)
       python 02_outcome_free_evaluation.py aggregate  (combine seeds and record SHA-256 hashes)
"""
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
import direct_vs_separate as D  # noqa: E402
from survival_tools import RESULTS, WORK, sha256  # noqa: E402

TAU = 365.25 * 5
BOOT = 999
SEEDS = (42, 52, 62, 72)
LABELS = {"ADOPT": "CANDIDATE SUPPORTED", "KEEP_REFERENCE": "REFERENCE SUPPORTED", "DEFER": "DEFER"}
INPUTS = WORK / "rotterdam_temporal" / "analysis_inputs"
OUT = RESULTS / "rotterdam_temporal"


def run(seed):
    d = INPUTS / f"seed{seed}"
    cal = pd.read_csv(d / "calibration.csv")
    target = pd.read_csv(d / "target.csv")
    if set(target.columns) != {"q0", "q1"}:
        raise RuntimeError("target inputs must contain only q0 and q1")
    rows, meta = D.analyze("Rotterdam-temporal", "candidate_vs_reference", seed, cal, target, TAU, BOOT,
                           rng_seed=D.stable_seed("temporal-rotterdam-20261110", seed, BOOT))
    OUT.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(rows)
    df["decision_direct"] = df.action_dir.map(LABELS); df["decision_separate"] = df.action_sep.map(LABELS)
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


if __name__ == "__main__":
    aggregate() if sys.argv[1] == "aggregate" else run(int(sys.argv[1]))
