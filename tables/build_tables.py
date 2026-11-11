"""Assemble the summary tables from the analysis outputs (see tables/README.md).

Requires results/ (distributed) and, for the calibration summaries of the applications, the stage-1
inputs under work/ (created by the application scripts after the data export).
"""
import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
RES = ROOT / "results"; WORK = ROOT / "work"; H = ROOT / "tables"
GAM = (0.0, 0.1, 0.2, 0.3, 0.5)
LAB = {"ADOPT": "CANDIDATE SUPPORTED", "KEEP_REFERENCE": "REFERENCE SUPPORTED", "DEFER": "DEFER",
       "CAND": "CANDIDATE SUPPORTED", "REF": "REFERENCE SUPPORTED"}
# Population % of source patients censored at or before tau without an event (Monte Carlo, n = 400,000).
BEFORE_TAU = {70: 64.9, 50: 42.5, 30: 23.5}


def wilson(k, n, z=1.959963984540054):
    p = k / n; d = 1 + z * z / n; c = (p + z * z / (2 * n)) / d; h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return c - h, c + h


def simulation_tables():
    cells, byg = [], []
    for nS in (100, 200, 300, 400, 500):
        for pct in (70, 50, 30):
            recs = [json.loads(l) for l in (RES / f"simulation/records/cell_{nS}_{pct}.jsonl").read_text().splitlines()]
            ok = [r for r in recs if not r["fatal"]]; n = len(ok)
            Wd = np.array([r["dir_width"] for r in ok]); Ws = np.array([r["sep_width"] for r in ok])
            row = {"n_source": nS, "censoring_target_pct": pct, "censored_before_tau_pct_population": BEFORE_TAU[pct],
                   "censoring_realized_mean": round(float(np.mean([r["cens"] for r in ok])), 4),
                   "events_mean": round(float(np.mean([r["events"] for r in ok])), 1), "replications": len(recs), "fatal": len(recs) - n}
            for m in ("dir", "sep"):
                k = sum(r[f"{m}_outer"] for r in ok); lo, hi = wilson(k, n)
                row.update({f"{m}_covered_n": k, f"{m}_coverage": round(k / n, 3), f"{m}_cov_wilson_low": round(lo, 3), f"{m}_cov_wilson_high": round(hi, 3)})
            row["id_ratio_mean_of_ratios_g>0"] = round(float((Ws[:, 1:] / Wd[:, 1:]).mean()), 3)
            row["id_ratio_ratio_of_means_g>0"] = round(float(Ws[:, 1:].mean() / Wd[:, 1:].mean()), 3)
            cells.append(row)
            for j, g in enumerate(GAM):
                b = {"n_source": nS, "censoring_target_pct": pct, "gamma": g, "replications": n,
                     "dir_width_mean": float(Wd[:, j].mean()), "sep_width_mean": float(Ws[:, j].mean()),
                     "id_ratio_mean_of_ratios": float((Ws[:, j] / Wd[:, j]).mean()) if g > 0 else None,
                     "id_ratio_ratio_of_means": float(Ws[:, j].mean() / Wd[:, j].mean()) if g > 0 else None}
                for m in ("dir", "sep"):
                    dec = [r[f"{m}_decision"][j] for r in ok]
                    b.update({f"{m}_correct_REF_n": dec.count("REF"), f"{m}_wrong_CAND_n": dec.count("CAND"), f"{m}_DEFER_n": dec.count("DEFER")})
                byg.append(b)
    pd.DataFrame(cells).to_csv(H / "simulation_conditions.csv", index=False)
    pd.DataFrame(byg).to_csv(H / "simulation_by_gamma.csv", index=False)


def application_tables():
    fx, cfg = [], []

    def add_rows(app, s, rows):
        for r in rows:
            fx.append({"application": app, "configuration": s, "gamma": r["gamma"], "L_dir": r["L_dir"], "U_dir": r["U_dir"],
                       "L_sep": r["L_sep"], "U_sep": r["U_sep"], "env_dir_lower": r["dir_sim_lower"], "env_dir_upper": r["dir_sim_upper"],
                       "env_sep_lower": r["sep_sim_lower"], "env_sep_upper": r["sep_sim_upper"],
                       "id_width_ratio": r["ID_width_ratio"] if np.isfinite(r["ID_width_ratio"]) else None,
                       "envelope_width_ratio": r["band_ratio"],
                       "decision_direct": LAB[r["action_dir"]], "decision_separate": LAB[r["action_sep"]]})

    def cal_summary(app, s, cal, tq0, tq1, tau, retro):
        n_target = len(tq0)
        same_side = round(100 * float((((np.asarray(tq0) - 0.5) * (np.asarray(tq1) - 0.5)) > 0).mean()), 1)
        cfg.append({"application": app, "configuration": s, "n_cal": len(cal), "cal_events_total": int(cal.event.sum()),
                    "cal_events_within_tau": int(((cal.time <= tau) & (cal.event == 1)).sum()),
                    "cal_censored_before_tau_pct": round(100 * float(((cal.time <= tau) & (cal.event == 0)).mean()), 1),
                    "n_target": n_target, "same_side_pct": same_side, "retro_marginalKM_contrast": retro.loc[s, "contrast"],
                    "retro_marginalKM_low": retro.loc[s, "ci_low"], "retro_marginalKM_high": retro.loc[s, "ci_high"]})

    pub = pd.read_csv(RES / "rotterdam_gbsg/direct_vs_separate.csv")
    rg = pd.read_csv(RES / "rotterdam_gbsg/retrospective_target_contrast.csv").set_index("seed")
    for s in (42, 52, 62, 72):
        add_rows("Rotterdam-GBSG", s, pub[pub.seed == s].sort_values("gamma").to_dict("records"))
        base = WORK / f"rotterdam_gbsg/analysis_inputs/seed{s}"
        cal = pd.read_csv(base / "calibration_GBSG.csv"); tg = pd.read_csv(base / "target_GBSG.csv")
        cal_summary("Rotterdam-GBSG", s, cal, tg.q0_surv5y, tg.q1_surv5y, 365.25 * 5, rg)
    for app, folder, tau in (("Rotterdam-temporal", "rotterdam_temporal", 365.25 * 5), ("NWTS3-NWTS4", "nwts", 365.25 * 3)):
        rt = pd.read_csv(RES / folder / "retrospective_target_contrast.csv").set_index("seed")
        for s in (42, 52, 62, 72):
            add_rows(app, s, pd.read_csv(RES / folder / f"seed{s}.csv").sort_values("gamma").to_dict("records"))
            cal = pd.read_csv(WORK / folder / f"analysis_inputs/seed{s}/calibration.csv")
            tg = pd.read_csv(WORK / folder / f"analysis_inputs/seed{s}/target.csv")
            cal_summary(app, s, cal, tg.q0, tg.q1, tau, rt)
    pd.DataFrame(fx).to_csv(H / "applications_fixed_configurations.csv", index=False)
    pd.DataFrame(cfg).to_csv(H / "applications_configuration_summary.csv", index=False)


def repeated_split_tables():
    rs_rows, rs_retro = [], []
    for app in ("gbsg", "temporal", "nwts"):
        recs = [json.loads(l) for l in (RES / f"repeated_splits/records/{app}.jsonl").read_text().splitlines()]
        ok = [r for r in recs if not r["fit_failure"]]
        retro = {json.loads(l)["split"]: json.loads(l) for l in (RES / f"repeated_splits/records/{app}_retro.jsonl").read_text().splitlines()}
        ev = [r["cal_events_tau"] for r in ok]
        for j, g in enumerate(GAM):
            C = [(r["split"], r["cells"][j]) for r in ok]
            row = {"application": app, "gamma": g, "splits": len(recs), "valid_splits": len(ok), "fit_failures": len(recs) - len(ok),
                   "cal_events_within_tau_median": float(np.median(ev)), "cal_events_within_tau_min": int(min(ev)), "cal_events_within_tau_max": int(max(ev)),
                   "dir_support_n": sum(c["action_dir"] != "DEFER" for _, c in C), "sep_support_n": sum(c["action_sep"] != "DEFER" for _, c in C),
                   "dir_only_support_n": sum(c["action_dir"] != "DEFER" and c["action_sep"] == "DEFER" for _, c in C),
                   "sep_only_support_n": sum(c["action_sep"] != "DEFER" and c["action_dir"] == "DEFER" for _, c in C),
                   "dir_REF_n": sum(c["action_dir"] == "REF" for _, c in C), "dir_CAND_n": sum(c["action_dir"] == "CAND" for _, c in C)}
            if g > 0:
                idr = np.array([c["ID_width_ratio"] for _, c in C]); bdr = np.array([c["band_ratio"] for _, c in C])
                row.update({"id_ratio_median": float(np.median(idr)), "id_ratio_q1": float(np.quantile(idr, .25)), "id_ratio_q3": float(np.quantile(idr, .75)),
                            "id_ratio_min": float(idr.min()), "id_ratio_max": float(idr.max()),
                            "env_ratio_median": float(np.median(bdr)), "env_ratio_q1": float(np.quantile(bdr, .25)), "env_ratio_q3": float(np.quantile(bdr, .75)),
                            "env_ratio_min": float(bdr.min()), "env_ratio_max": float(bdr.max())})
            for m in ("dir", "sep"):
                cert = [(s, c[f"action_{m}"]) for s, c in C if c[f"action_{m}"] != "DEFER"]; inc = opp = 0
                for s, a in cert:
                    rr = retro[s]; sgn = -1 if a == "CAND" else 1
                    if np.sign(rr["contrast"]) != sgn:
                        if (sgn == -1 and rr["ci_low"] > 0) or (sgn == 1 and rr["ci_high"] < 0): opp += 1
                        else: inc += 1
                row.update({f"{m}_disagree_retro_ci_includes0_n": inc, f"{m}_disagree_retro_ci_opposite_n": opp})
            rs_rows.append(row)
        rr = pd.DataFrame(retro.values())
        rs_retro.append({"application": app, "splits": len(rr), "contrast_median": rr.contrast.median(), "contrast_min": rr.contrast.min(),
                         "contrast_max": rr.contrast.max(), "ci_entirely_below0_n": int((rr.ci_high < 0).sum()),
                         "ci_includes0_n": int(((rr.ci_low <= 0) & (rr.ci_high >= 0)).sum()), "ci_entirely_above0_n": int((rr.ci_low > 0).sum())})
    pd.DataFrame(rs_rows).to_csv(H / "repeated_splits_by_gamma.csv", index=False)
    pd.DataFrame(rs_retro).to_csv(H / "repeated_splits_retrospective.csv", index=False)


if __name__ == "__main__":
    which = sys.argv[1:] or ["simulation", "applications", "repeated_splits"]
    if "simulation" in which: simulation_tables()
    if "applications" in which: application_tables()
    if "repeated_splits" in which: repeated_split_tables()
    print("tables written to", H)
