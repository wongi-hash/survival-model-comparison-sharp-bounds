# Simulation factorial v2 (identification, inference, decision) — fixed 2026-11-10 before any run

Purpose: replace the application-matched factorial (n_S 114/516, censoring 65/50/38%, n_T 80) by a
generic grid and evaluate, on the same bootstrap draws, direct and separate identification at the
three levels used in the manuscript: identification width, simultaneous coverage, and decisions.

- Data-generating mechanism, estimator (conditional Cox-Breslow IPCW + IPCW-weighted logistic),
  joint pairs bootstrap and sup-t critical value: unchanged from `RUN_COVERAGE_FACTORIAL_nS_CENSORING.py`
  (true shift 0.2*tanh(x); target X ~ N(0.5, 1); reference truly superior at every gamma).
- Grid: n_S in {100, 200, 300, 400, 500} x source censoring in {70%, 50%, 30%}; n_T = 500;
  gamma in {0, 0.1, 0.2, 0.3, 0.5}; 1,000 replications per cell; B = 199.
- Seeds: replication seed stable_seed("factorial-v2", n_S, round(100*censoring), rep); bootstrap
  seed = replication seed + 800,000,003. Direct and separate use the same draws; each method uses its
  own sup-t critical value.
- Recorded per replication: realized censoring and events; direct and separate identified-set widths
  at each gamma; outer coverage of each method's true population endpoints; decision at each gamma
  for each method (REFERENCE SUPPORTED = correct direction here; CANDIDATE SUPPORTED = wrong; DEFER).
- Reported for all 15 cells: mean width ratio, outer coverage (Wilson 95% CI), correct, wrong and
  DEFER rates by gamma. No sample-size threshold will be claimed.
