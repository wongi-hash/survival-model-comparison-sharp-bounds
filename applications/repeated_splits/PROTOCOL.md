# Repeated-split sensitivity analysis — protocol fixed 2026-11-10 before any run

Retrospective, pre-planned split-sensitivity analysis. Target outcomes of all three applications had
already been examined for the four primary configurations; this analysis does not select splits.

## Aim
Quantify how the width gain and the decision differences of direct versus separate identification
depend on the development/calibration split of the source data.

## Applications (definitions identical to the primary analyses)
- GBSG: source node-positive Rotterdam (all years, n = 1,546), target GBSG (n = 686), tau = 5 years;
  reference = age^3, age^3 log(age), menopause, size, hormonal therapy; candidate adds nodes^(-1/2), ER.
- Temporal: source node-positive Rotterdam operated <= 1988 (n = 860), target 1989-1993 (n = 686),
  tau = 5 years; same predictors.
- NWTS: source NWTS-3 (n = 1,857), target NWTS-4 (n = 2,171), tau = 3 years; reference = stage,
  age, institutional histology; candidate adds central histology.

## Splits and implementation
- 100 splits per application with split seeds 1001-1100, fixed in advance; event-stratified 2:1
  development/calibration split (numpy). The four primary configurations are not part of this set.
- Unified predictor implementation (`cox_unified.py`: ridge Cox, Efron ties, interpolated Breslow
  baseline), which reproduces the published Rotterdam-GBSG predictions to ~1e-5.
- Outcome-free evaluation: `direct_vs_separate.analyze`, B = 999 (as in the primary analyses),
  bootstrap seed stable_seed("repsplit-20261110", application, split, 999).

## Recorded and reported
- Per split and gamma in {0, 0.1, 0.2, 0.3, 0.5}: identified-set and envelope widths and their ratios;
  decisions of both methods; model-fitting failures recorded separately from DEFER.
- Summaries across splits: distribution (median, IQR, range) of width ratios; proportion of splits
  in which each method certifies a direction; proportion in which only the direct method certifies.
- Retrospective stage (run only after the outcome-free results are saved and hashed): per split, the
  target IPCW Brier contrast at tau with a 999-resample percentile interval. Reported as the
  proportion of certified cells whose direction disagrees with the retrospective point estimate,
  separately for (a) retrospective interval containing zero and (b) retrospective interval entirely on
  the opposite side of zero. Not described as an error rate.
- Splits reuse the same patients; no binomial confidence intervals are attached to proportions
  across splits, and the splits are not described as independent validations.
