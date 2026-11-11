"""Stage 2: outcome-free direct versus separate evaluation for the Rotterdam -> GBSG application.

Uses the archived bootstrap seeds in src/archived_bootstrap_seeds.json, so a rerun reproduces the
reported bootstrap draws. Writes results/rotterdam_gbsg/direct_vs_separate.csv.
"""
from __future__ import annotations
import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
import direct_vs_separate as core  # noqa: E402

INPUTS = ROOT / "work" / "rotterdam_gbsg" / "analysis_inputs"
OUT = ROOT / "results" / "rotterdam_gbsg"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bootstrap", type=int, default=999)
    args = parser.parse_args()
    rows, metadata = [], []
    for seed in core.SEEDS:
        folder = INPUTS / f"seed{seed}"
        calibration = pd.read_csv(folder / "calibration_GBSG.csv").rename(columns={"q0_surv5y": "q0", "q1_surv5y": "q1"})
        target = pd.read_csv(folder / "target_GBSG.csv").rename(columns={"q0_surv5y": "q0", "q1_surv5y": "q1"})
        if {"time", "event"} & set(target):
            raise RuntimeError("target inputs must not contain outcomes")
        result, meta = core.analyze("Rotterdam-GBSG", "supporting_clinical_candidate", seed,
                                    calibration, target, 365.25 * 5, args.bootstrap)
        rows.extend(result)
        metadata.append(meta)
    OUT.mkdir(parents=True, exist_ok=True)
    output = OUT / "direct_vs_separate.csv"
    pd.DataFrame(rows).to_csv(output, index=False)
    manifest = {"status": "PASS", "created_utc": datetime.now(timezone.utc).isoformat(),
                "target_outcomes_read": False, "bootstrap_draws": args.bootstrap,
                "core_sha256": core.sha256(ROOT / "src" / "direct_vs_separate.py"),
                "output_sha256": core.sha256(output), "jobs": metadata}
    (OUT / "OUTCOME_FREE_MANIFEST.json").write_text(json.dumps(manifest, indent=2, default=str), encoding="utf-8")
    print("[PASS]", len(rows), "rows ->", output)


if __name__ == "__main__":
    main()
