# Cross-study application: NWTS-3 -> NWTS-4 (Wilms tumour)
Protocol fixed on 2026-11-10 before any target (NWTS-4) outcome was examined.

## Question
Two locked relapse-prediction models are developed on NWTS-3 patients. Before three-year relapse
outcomes are available for patients treated in the next study (NWTS-4), does adding central-pathology
histology to a model using institutional histology improve prediction of three-year relapse-free
survival in NWTS-4?

## Data
- `survival::nwtco` exported from R to `nwtco.csv` (SHA-256 b5fa39f1...7549), n = 4,028; full cohort
  (the case-cohort indicator is not used).
- Source: study = 3 (NWTS-3, n = 1,857). Target: study = 4 (NWTS-4, n = 2,171).
- Inspected before this protocol: covariate distributions by study; source (NWTS-3) relapse counts by
  year (267 relapses and 72 censorings before 3 years). Target outcomes were not examined.

## Outcome and horizon
- Relapse (rel) with time edrel (days); administratively truncated at 5 years; tau = 3 years
  (1,095.75 days).

## Locked predictors
- Seeds 42, 52, 62, 72; event-stratified 2:1 split of the source into predictor development and
  non-overlapping source calibration samples (same numpy routine as the temporal Rotterdam application).
- Reference q0: Cox model on stage II, III, IV indicators (vs I), age in years, institutional
  unfavourable histology (instit = 2).
- Candidate q1: q0 covariates plus central-pathology unfavourable histology (histol = 2).
- Ridge Cox (penalty 0.01, standardized covariates, Breslow baseline); q = predicted S(tau | x).

## Outcome-free evaluation (unchanged primary procedure)
- `code/inference/direct_vs_separate.py::analyze`, B = 999, sup-t over gamma in {0, 0.01, ..., 1};
  reported gamma in {0, 0.1, 0.2, 0.3, 0.5}; decisions CANDIDATE SUPPORTED / REFERENCE SUPPORTED / DEFER.
- Target inputs contain only q0 and q1. All 20 cells reported regardless of result.

## Retrospective comparison (after outcome-free results are saved and hashed)
- Target IPCW Brier contrast at tau with marginal Kaplan-Meier censoring; 999 bootstrap percentile
  interval; smallest gamma whose direct identified set contains the retrospective contrast.
  Descriptive only.

## Amendment 2 (2026-11-10, implementation unification)
The predictor implementation was changed from Breslow ties with a step-function baseline at tau to
Efron ties with a linearly interpolated cumulative baseline hazard, the convention of the lifelines
pipeline used for Rotterdam-GBSG. With identical splits this reproduces the published Rotterdam-GBSG
predictions to about 1e-5, IPCW weights to about 1e-6, and identified-set endpoints to about 3e-6.
All stages were rerun; results from the first implementation are archived in `*_v0_breslow_step/`.
