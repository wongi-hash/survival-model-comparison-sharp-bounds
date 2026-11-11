# Verification of the released code

The released scripts were re-executed on 10 November 2026 in a second computing environment (Python 3.10, NumPy 2.2.6,
pandas 2.3.3, Matplotlib 3.10.9) with the data exported as in `data/README.md` (all three SHA-256 hashes
matched), and compared with the distributed results.

| Analysis | Re-executed | Comparison with distributed results |
|---|---|---|
| Rotterdam earlier → later | Stages 1–3, all four configurations (B = 999) | Stage-1 input hashes identical; all result files byte-identical |
| NWTS-3 → NWTS-4 | Stages 1–3, all four configurations (B = 999) | Outcome-free result files byte-identical; retrospective contrasts identical |
| Rotterdam → GBSG | Stage 1 (separately, see below); stages 2–3 from the archived stage-1 inputs (B = 999) | Stage-1 predictions identical; stage-2 values within 9e-16 with identical decisions; retrospective contrasts identical |
| Repeated splits | One split per application (1001, 1050, 1100): outcome-free record and retrospective point estimate | Outcome-free records identical for two applications and within 8e-15 for the third; decisions identical; retrospective point estimates within 3e-17 |
| Simulation | Replications 0–2 of two conditions and 0–1 of a third | Records identical or within 1.1e-14; decisions and coverage indicators identical |
| Tables | `tables/build_tables.py` | All values identical to those reported; the application tables add the same-side percentage and leave the undefined width ratio at gamma = 0 blank |
| Figure 2 | `figures/make_figure2.py` | Regenerated from the distributed results |

Differences of order 1e-15 arise from floating-point arithmetic on different computers.

Stage 1 of the Rotterdam → GBSG application (predictors fitted with lifelines) was re-executed separately
on Windows with Python 3.14 and the packages in `requirements.txt`. The predicted survival probabilities in
the calibration and target inputs were identical (maximum absolute difference 0) to the archived stage-1
inputs, which had been produced under Python 3.8.5 and were used for the stage 2–3 checks above.
Recomputing this application with the numpy predictor implementation in `src/cox_ridge.py` changed
envelope limits by at most 2.2e-6, with identical decisions.

The record files in `results/simulation/records/` and `results/repeated_splits/records/` and the
outcome-free result files of the two protocol-based applications match the hashes in
`provenance/hash_records/`.
