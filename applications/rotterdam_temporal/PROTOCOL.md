# Temporal deployment application: Rotterdam early -> Rotterdam late
Protocol fixed on 2026-11-10 before any target outcome was examined.

## Question
A breast-cancer centre has two locked prognostic models developed on patients operated in
1978-1988. Before five-year outcomes of patients operated in 1989-1993 have matured, does the
candidate model (reference + estrogen-receptor status) predict five-year recurrence-free survival
more accurately than the reference model in the later cohort?

## Data
- `survival::rotterdam` (public; same file as the Rotterdam-GBSG pipeline), node-positive patients
  (nodes > 0), n = 1,546.
- Source: year of surgery <= 1988 (n = 860). Target: year of surgery 1989-1993 (n = 686).
- Only covariates were inspected before this protocol (counts by year; source event count 560 by
  seven years; hormonal therapy use 7% vs 33% before/after 1987). Target outcomes were not examined.

## Outcome and horizon
- Recurrence-free survival (recurrence or death), administratively truncated at 7 years; tau = 5 years.
  Definitions identical to `code/rotterdam_gbsg/common.py::source_outcome`.

## Locked predictors
- Seeds 42, 52, 62, 72. Within the source, event-stratified 2:1 split into predictor development and
  non-overlapping source calibration samples.
- Reference q0: Cox model on age^3, age^3*log(age), menopausal status, tumour size (20-50, >50 mm),
  hormonal therapy. Candidate q1: same plus nodes^(-1/2) and estrogen receptor.
  (Columns identical to `01_prepare_locked_predictors.py` of the Rotterdam-GBSG pipeline.)
  Amendment 1 (2026-11-10, before any analysis was run): the first draft listed nodes^(-1/2) in the
  reference model; corrected to match the published pipeline.
- Ridge Cox partial likelihood (penalty 0.01 on standardized covariates, mean log-likelihood scale,
  Breslow baseline), implemented in numpy because lifelines is unavailable in the analysis environment.
  q = predicted S(tau | x).

## Outcome-free evaluation (unchanged primary procedure)
- `code/inference/direct_vs_separate.py::analyze`: conditional Cox-Breslow IPCW given Z = (q0, q1),
  IPCW-weighted logistic outcome model, joint pairs bootstrap B = 999, sup-t critical value over
  gamma in {0, 0.01, ..., 1}; reported gamma in {0, 0.1, 0.2, 0.3, 0.5}.
- Decisions: CANDIDATE SUPPORTED if the direct (or separate) simultaneous envelope lies below 0,
  REFERENCE SUPPORTED if above 0, otherwise DEFER.
- Target data passed to the procedure contain only q0 and q1.
- All 20 cells (4 seeds x 5 gamma) are reported regardless of result.

## Retrospective comparison (only after the outcome-free results are saved and hashed)
- Target Brier-risk contrast at tau estimated by IPCW with a marginal Kaplan-Meier censoring
  distribution in the target; 999 bootstrap resamples, percentile 95% interval.
- Purpose: descriptive check of whether certified directions agree with the later-observed contrast.
  It does not alter any outcome-free decision.

## Amendment 2 (2026-11-10, implementation unification)
The predictor implementation was changed from Breslow ties with a step-function baseline at tau to
Efron ties with a linearly interpolated cumulative baseline hazard, the convention of the lifelines
pipeline used for Rotterdam-GBSG. With identical splits this reproduces the published Rotterdam-GBSG
predictions to about 1e-5, IPCW weights to about 1e-6, and identified-set endpoints to about 3e-6.
All stages were rerun; results from the first implementation are archived in `*_v0_breslow_step/`.
