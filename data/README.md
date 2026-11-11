# Data

No patient-level data are distributed with this repository. The three datasets are public and are
distributed with the R package `survival` (Therneau TM. A package for survival analysis in R.
R package version 3.8-3, 2024. https://CRAN.R-project.org/package=survival).

| File | R object | Use |
|---|---|---|
| `rotterdam.csv` | `survival::rotterdam` | Rotterdam → GBSG (source); Rotterdam earlier → later (source and target) |
| `gbsg.csv` | `survival::gbsg` | Rotterdam → GBSG (target) |
| `nwtco.csv` | `survival::nwtco` | NWTS-3 → NWTS-4 (source: study 3; target: study 4) |

## Export

```bash
cd data
Rscript export_survival_data.R
python prepare_data.py
```

`prepare_data.py` normalizes line endings and checks the SHA-256 hashes of the files used in the
analyses:

| File | SHA-256 |
|---|---|
| `rotterdam.csv` | `62703670d3be5d49c5476fb7131b020bdca9f1c25bb4522fd6abe8f9dab61f1b` |
| `gbsg.csv` | `9fa0fef0575d04d3273869be28405307f1b02993751498e8a0c844541cccd5d9` |
| `nwtco.csv` | `b5fa39f11d2adb8f51845dd302806da7ed09c4eecb453548a70f0692ebf77549` |

The analysis scripts stop if a hash differs. A difference may arise from a different version of the
`survival` package or of the CSV writer; in that case compare the exported values with the package
documentation before changing any expected hash.
