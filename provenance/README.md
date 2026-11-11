# Provenance

## Hash records (`hash_records/`)

Files written during the analyses. They list SHA-256 hashes with UTC time stamps. Except for the
Rotterdam → GBSG record (see below), they are reproduced unchanged.

| File | Content |
|---|---|
| `*_PROTOCOL.sha256` | Hashes of each protocol (and, where listed, runner scripts) at the time they were fixed or amended |
| `rotterdam_temporal_OUTCOME_FREE_SHA256.txt`, `nwts_OUTCOME_FREE_SHA256.txt` | Outcome-free results, hashed before the retrospective stage |
| `repeated_splits_OUTCOME_FREE_SHA256.txt` | Outcome-free split records, hashed before the retrospective stage |
| `repeated_splits_RETRO_SHA256.txt` | Retrospective split records |
| `simulation_RESULTS_SHA256.txt` | Simulation records |
| `rotterdam_gbsg_stage1_2026-09-02.json` | Hashes of the Rotterdam → GBSG stage-1 inputs recorded on 2 September 2026; fields extracted from the original stage-1 manifest with local file paths removed |

The simulation and repeated-split record files distributed in `results/`, and the outcome-free result files
of the Rotterdam earlier → later and NWTS-3 → NWTS-4 applications, match these hashes. File names and script
names in the records refer to the layout used during the analyses.

## Timeline (UTC)

| Date and time | Event |
|---|---|
| 2 Sep 2026, 05:21 | Rotterdam → GBSG predictor definitions and settings fixed; stage-1 inputs created and hashed |
| 10 Nov 2026, 01:36 | Rotterdam earlier → later protocol fixed; Amendment 1 at 01:36 |
| 10 Nov 2026, 01:38 | Rotterdam earlier → later initial outcome-free outputs hashed; target outcomes then read |
| 10 Nov 2026, 01:49 | NWTS-3 → NWTS-4 protocol fixed |
| 10 Nov 2026, 01:52 | NWTS-3 → NWTS-4 initial outcome-free outputs hashed; target outcomes then read |
| 10 Nov 2026, 01:59 | Simulation protocol fixed |
| 10 Nov 2026, 02:15–02:16 | Amendment 2 (predictor implementation) for both applications |
| 10 Nov 2026, 02:16 and 02:18 | Rerun outcome-free outputs hashed (Rotterdam earlier → later; NWTS-3 → NWTS-4), followed by the retrospective stage |
| 10 Nov 2026, 02:21 | Repeated-split protocol fixed |
| 10 Nov 2026, 03:51 | Repeated-split outcome-free records hashed; retrospective records hashed at 03:58 |
| 10 Nov 2026, 04:40 | Simulation records (15 conditions x 1,000 replications) hashed |

The analyses after 2 September 2026 were run within one day with fixed seeds; the simulation replications
were run in parallel on two computers.

## Analysis history

**Rotterdam earlier → later and NWTS-3 → NWTS-4.** Written protocols were fixed before any target outcome
was examined; covariate distributions and source event counts had been inspected beforehand. Amendment 1
(Rotterdam earlier → later), made before any analysis was run, corrected the reference covariate list.
The initial analysis used a predictor implementation with Breslow ties and a step-function baseline at
tau; its outcome-free outputs were saved and hashed, and target outcomes were then read for the
retrospective comparison. Amendment 2, made after that first access to target outcomes, replaced the
predictor implementation by the one in `src/cox_ridge.py` (Efron ties and a linearly interpolated baseline
cumulative hazard at tau), which reproduces the Rotterdam → GBSG predictors. All stages were rerun, and the
new outcome-free outputs were hashed before the retrospective stage was rerun. Compared with the initial
run, the directional decisions for both methods were identical in all 40 configuration-by-sensitivity
cells; identified-interval endpoints differed by at most 2.2e-4, envelope limits by at most 2.6e-4, and
retrospective contrasts by at most 4.0e-5. The initial results are in `initial_implementation/`.

**Rotterdam → GBSG.** The final predictor definitions and analysis settings were fixed on 2 September 2026,
before the finalized outcome-free computation. GBSG outcomes were excluded from that computation, although
they had previously been examined for a different predictor contrast.
`applications/rotterdam_gbsg/SPECIFICATION.md` restates these definitions and settings; no time-stamped hash
of a specification document from September 2026 is distributed. The time-stamped record of that date is
`hash_records/rotterdam_gbsg_stage1_2026-09-02.json`, whose hashes match the stage-1 inputs used for the
reported results. A rerun of stage 1 in November 2026 (Windows, Python 3.14) reproduced the predicted
survival probabilities exactly (see `VERIFICATION.md`).

**Repeated splits.** The protocol was fixed before the analysis was run, after the target outcomes of the
four primary configurations had been examined; it does not select splits. Outcome-free records of all 300
splits were hashed before the retrospective stage.

**Simulation.** The protocol was fixed before the run. Its description of the data-generating mechanism
mentions a true shift of 0.2·tanh(x), carried over from the description of the earlier application-matched
design (n_S = 114/516, n_T = 80) that this factorial replaced; the shift does not enter the code. Because
coverage and decision correctness are defined against population identified-interval endpoints, target
outcomes were not generated. Records were produced on two computers; replications recomputed on a single
computer agree to within 1.1e-14 with identical decisions and coverage indicators (see `VERIFICATION.md`).
