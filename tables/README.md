# Summary tables

Built by `build_tables.py` from `results/` (and, for calibration summaries, the stage-1 inputs under `work/`).

| File | Content |
|---|---|
| `simulation_conditions.csv` | 15 conditions: source size, target and realized overall censoring, population proportion censored before tau, mean source events, simultaneous coverage with Wilson 95% interval (direct, separate), identification-width ratio |
| `simulation_by_gamma.csv` | 15 conditions x 5 sensitivity values: mean identified-interval widths, width ratios, counts of REFERENCE SUPPORTED (correct), CANDIDATE SUPPORTED (incorrect) and DEFER for each method (1,000 replications) |
| `applications_fixed_configurations.csv` | 3 applications x 4 configurations x 5 sensitivity values: direct and separate identified intervals, simultaneous envelopes, width ratios, decisions |
| `applications_configuration_summary.csv` | Per configuration: calibration size, calibration events (total and by tau), % censored before tau, target size, % of target profiles with both predictions on the same side of 0.5, retrospective contrast with 95% interval |
| `repeated_splits_by_gamma.csv` | 3 applications x 5 sensitivity values: valid splits, fitting failures, calibration events by tau, support counts (direct, separate, direct only, separate only), width-ratio median/IQR/range, disagreement counts |
| `repeated_splits_retrospective.csv` | Per application: retrospective contrast median and range across splits; numbers of intervals below, including and above zero |

## Definitions

- **Decision labels.** In the per-configuration outputs under `results/`, `ADOPT` and `KEEP_REFERENCE`
  (and `CAND` and `REF` in the record files) denote CANDIDATE SUPPORTED and REFERENCE SUPPORTED.

- **Identification-width ratio (simulation).** `id_ratio_mean_of_ratios`: separate width divided by direct
  width in each replication, averaged over replications (and over gamma > 0 in `simulation_conditions.csv`).
  `id_ratio_ratio_of_means`: mean separate width divided by mean direct width.
- **Width ratios (applications, repeated splits).** Separate width divided by direct width for the
  identified interval or the simultaneous envelope. The identified-interval ratio is undefined at gamma = 0,
  where both intervals are points, and is left blank in `applications_fixed_configurations.csv` (the record
  files store it as Infinity).
- **Coverage (simulation).** The method's 95% simultaneous envelope contains its own population
  identified-interval endpoints at every gamma in {0, 0.1, 0.2, 0.3, 0.5}.
- **Correct and incorrect support (simulation).** The reference predictor has the lower target Brier risk
  for every target outcome distribution in the sensitivity class with gamma <= 0.5, so REFERENCE SUPPORTED
  is correct and CANDIDATE SUPPORTED is incorrect.
- **Disagreement (repeated splits).** Among splits in which a method supported a direction at a given
  gamma, the number whose direction was opposite to the sign of the retrospective point estimate, split by
  whether the retrospective interval included zero or lay on the opposite side.
- **Events.** Applications: `cal_events_within_tau` counts events at or before tau in the calibration sample;
  `cal_events_total` counts events over follow-up (administratively censored at seven years for Rotterdam
  and five years for NWTS). Simulation: `events_mean` is the mean number of observed source events.
- **Censored before tau.** Applications: proportion of calibration patients censored at or before tau
  without an event. Simulation: population value (Monte Carlo, n = 400,000).
- **Retrospective contrast.** Candidate-minus-reference target Brier risk at tau estimated with IPCW and a
  marginal Kaplan-Meier censoring distribution; percentile 95% interval from 999 resamples of target
  patients. Computed after the outcome-free analysis and not used by the comparison procedure.
