"""Stage 3 (run only after stage 2 outputs are saved and hashed): retrospective target Brier-risk
difference at tau, IPCW with a marginal Kaplan-Meier censoring distribution, 999-resample percentile
interval. Descriptive only; it does not change any outcome-free decision.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from survival_tools import RESULTS, WORK, ipcw_brier_scores, retrospective_interval  # noqa: E402

TAU = 365.25 * 5
BOOT = 999
SEEDS = (42, 52, 62, 72)
W = WORK / "rotterdam_temporal"
OUT = RESULTS / "rotterdam_temporal"


def main():
    if not (OUT / "OUTCOME_FREE_SHA256.txt").is_file():
        raise RuntimeError("run '02_outcome_free_evaluation.py aggregate' first")
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
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "retrospective_target_contrast.csv", index=False)
    print(df.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
