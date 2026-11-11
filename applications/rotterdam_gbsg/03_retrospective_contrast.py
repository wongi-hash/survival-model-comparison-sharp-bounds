"""Stage 3 (run only after stage 2): retrospective target Brier-risk difference at five years, IPCW with
a marginal Kaplan-Meier censoring distribution, 999-resample percentile interval. Descriptive only.
"""
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from survival_tools import ipcw_brier_scores, retrospective_interval  # noqa: E402

TAU = 365.25 * 5
BOOT = 999
SEEDS = (42, 52, 62, 72)
W = ROOT / "work" / "rotterdam_gbsg"
OUT = ROOT / "results" / "rotterdam_gbsg"


def main():
    if not (OUT / "direct_vs_separate.csv").is_file():
        raise RuntimeError("run 02_outcome_free_evaluation.py first")
    oc = pd.read_csv(W / "sealed_target_outcomes" / "target_GBSG_outcomes.csv")
    rows = []
    for seed in SEEDS:
        tg = pd.read_csv(W / "analysis_inputs" / f"seed{seed}" / "target_GBSG.csv")
        if not (tg.patient_id.astype(str).to_numpy() == oc.patient_id.astype(str).to_numpy()).all():
            raise RuntimeError("patient order mismatch")
        t, e = oc.time.to_numpy(float), oc.event.to_numpy(int)
        q0, q1 = tg.q0_surv5y.to_numpy(), tg.q1_surv5y.to_numpy()
        b0, b1 = ipcw_brier_scores(t, e, q0, q1, TAU)
        est, lo, hi = retrospective_interval(t, e, q0, q1, TAU, seed=20261110 + seed, boot=BOOT)
        rows.append({"seed": seed, "brier_ref": b0, "brier_cand": b1, "contrast": est, "ci_low": lo, "ci_high": hi,
                     "target_censored_before_tau": float(((t <= TAU) & (e == 0)).mean())})
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "retrospective_target_contrast.csv", index=False)
    print(df.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
