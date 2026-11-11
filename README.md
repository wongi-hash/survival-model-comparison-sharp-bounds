# Comparing survival prediction models with sharp bounds on Brier-risk differences before local outcomes are available

Code and aggregate results for the simulation study and the three public-data applications in the
manuscript of the same title (Choi WG, Kim S-H, Kim MS).

Two fixed survival predictors, a reference `q0` and a candidate `q1`, are compared at a fixed horizon in a
target population whose outcomes are not yet available. Under a bounded log-odds shift `|delta(z)| <= gamma`
between source and target conditional survival given the prediction profile `Z = (q0, q1)`, the code computes

- sharp bounds for the candidate-minus-reference Brier-risk difference (direct identification);
- the interval obtained by bounding each risk separately and subtracting (separate-risk subtraction);
- censoring-adjusted estimates and joint-bootstrap simultaneous 95% confidence envelopes;
- the decisions CANDIDATE SUPPORTED (upper limit below zero), REFERENCE SUPPORTED (lower limit above zero)
  and DEFER.

## Contents

| Folder | Content |
|---|---|
| `src/` | Core estimator (`direct_vs_separate.py`), ridge Cox model, shared helpers, archived bootstrap seeds |
| `simulation/` | Data-generating mechanism, simulation runner, protocol |
| `applications/rotterdam_gbsg/` | Rotterdam → German Breast Cancer Study Group (other hospital), 5-year horizon |
| `applications/rotterdam_temporal/` | Rotterdam, surgery 1978–1988 → 1989–1993 (later patients), 5-year horizon |
| `applications/nwts/` | National Wilms Tumor Study 3 → 4 (subsequent study), 3-year horizon |
| `applications/repeated_splits/` | 100 development–calibration splits per application |
| `results/` | Aggregate outputs of all analyses (simulation records, per-configuration results, split records) |
| `tables/` | Summary tables built from `results/` (definitions in `tables/README.md`) |
| `figures/` | Figure 2 and its script |
| `data/` | Export script for the public datasets and expected SHA-256 hashes (no data distributed) |
| `provenance/` | Hash records made during the analyses, initial-implementation results, and analysis history |

Directory `work/` is created locally for patient-level intermediate files and is not distributed.

## Requirements

Python 3.9 or later with the packages in `requirements.txt`. R with the `survival` package is needed only
to export the data. `lifelines` and `scikit-learn` are needed only for stage 1 of the Rotterdam → GBSG
application.

```bash
pip install -r requirements.txt
```

## Reproduction

Run from the repository root. Outputs are written to `results/`, `tables/` and `figures/`, overwriting
the distributed versions; compare with the distributed files (for example with `git diff`).

```bash
# 0. Data (see data/README.md)
cd data && Rscript export_survival_data.R && python prepare_data.py && cd ..

# 1. Simulation (15 conditions x 1,000 replications; resumable; long-running)
python simulation/run_factorial.py all
python simulation/run_factorial.py summarize
#    single condition, e.g.:  python simulation/run_factorial.py 100 70

# 2. Rotterdam -> GBSG
python applications/rotterdam_gbsg/01_prepare_predictors.py
python applications/rotterdam_gbsg/02_outcome_free_evaluation.py
python applications/rotterdam_gbsg/03_retrospective_contrast.py

# 3. Rotterdam earlier -> later
python applications/rotterdam_temporal/01_prepare_predictors.py
for s in 42 52 62 72; do python applications/rotterdam_temporal/02_outcome_free_evaluation.py $s; done
python applications/rotterdam_temporal/02_outcome_free_evaluation.py aggregate
python applications/rotterdam_temporal/03_retrospective_contrast.py

# 4. NWTS-3 -> NWTS-4
python applications/nwts/run_nwts.py prepare
for s in 42 52 62 72; do python applications/nwts/run_nwts.py evaluate $s; done
python applications/nwts/run_nwts.py aggregate
python applications/nwts/run_nwts.py retrospective

# 5. Repeated splits (resumable; long-running)
for a in gbsg temporal nwts; do python applications/repeated_splits/run_repeated_splits.py run $a; done
python applications/repeated_splits/run_repeated_splits.py hash
for a in gbsg temporal nwts; do python applications/repeated_splits/run_repeated_splits.py retro $a; done
python applications/repeated_splits/run_repeated_splits.py hash-retro

# 6. Tables and Figure 2 (step 6 needs the stage-1 inputs under work/ from steps 2-4)
python tables/build_tables.py
python figures/make_figure2.py
```

In each application, stage 2 reads target inputs that contain only `q0` and `q1`; target outcomes are
written to a separate folder in stage 1 and are read only by stage 3.

## Reproducibility notes

- All random numbers are generated from fixed seeds. The Rotterdam → GBSG bootstrap seeds are stored in
  `src/archived_bootstrap_seeds.json`; the other analyses derive their seeds from fixed labels.
- `VERIFICATION.md` reports a re-execution of the released code against the distributed results.
- The protocols in `applications/*/PROTOCOL.md` and `simulation/PROTOCOL.md` are the documents fixed before
  the analyses and are reproduced unchanged; they refer to the file names used at that time.
  `applications/rotterdam_gbsg/SPECIFICATION.md` restates the Rotterdam → GBSG settings (see `provenance/`). The scripts were reorganized for release (file paths,
  imports and comments only), so their file hashes differ from the script hashes listed in
  `provenance/hash_records/` and in the stage-1 manifests under `results/`.

## License

Code: MIT License (see `LICENSE`). The datasets are not redistributed; they are available in the R package
`survival`.

## Citation

See `CITATION.cff`.
