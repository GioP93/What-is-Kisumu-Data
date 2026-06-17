# Between-laboratory variability in the Kisumu *Anopheles gambiae* reference strain

Anonymised data and analysis code for the study of phenotypic susceptibility
variation in the *Anopheles gambiae* "Kisumu" insecticide-susceptible
reference strain across sixteen laboratories, tested against four insecticides
(permethrin, alpha-cypermethrin, DDT and pirimiphos-methyl) in WHO tube
bioassays and analysed with a Bayesian five-parameter logistic
concentration–response model.

This repository lets you reproduce the LC50 estimates, between/within-laboratory
variability metrics, dose-response model checks, significance groupings, and the
cross-insecticide Spearman rank correlation reported in the paper.

## Contents

```
.
├── data/                         Anonymised input data (one file per insecticide)
│   ├── permethrin.csv
│   ├── alpha_cypermethrin.csv
│   ├── ddt.csv
│   └── pirimiphos_methyl.csv
├── stan/
│   └── example_model.stan        <-- YOU MUST ADD THIS FILE (see below)
├── scripts/
│   ├── 1_variability_analysis.R  Fits all four insecticides; LC50s, plots, metrics
│   └── 2_spearman_correlation.R  Cross-insecticide Spearman (run after script 1)
├── outputs/                      Created on first run (empty in the archive)
├── run_all.R                     Runs script 1 then script 2
├── kisumu_variability.Rproj      Open this first so paths resolve
└── README.md
```

## Site anonymisation

Testing laboratories are identified only by a letter code (A–P). The real
laboratory names are not included in this archive. The code `K_star` in the
data files corresponds to the display label `K*` used in the paper's figures
and tables. The mapping between letter codes and the figures/tables in the
paper is consistent throughout.

## Requirements

- R (≥ 4.2 recommended)
- A working C++ toolchain for `rstan` (Rtools on Windows; Xcode CLT on macOS)
- The following R packages (the scripts install/load them via `pacman`):
  `pacman`, `here`, `rio`, `tidyverse`, `reshape2`, `rstan`, `loo`, `boot`,
  `egg`, `gridExtra`, `ggpubr`, `ggsci`, `viridisLite`, `patchwork`

If `pacman` is not installed: `install.packages("pacman")`.

## How to run

1. Unzip the archive.
2. Add your `example_model.stan` to the `stan/` folder.
3. Open `kisumu_variability.Rproj` in RStudio (this sets the working
   directory so all paths resolve via the `here` package).
4. Either run everything:
   ```r
   source("run_all.R")
   ```
   or run the scripts individually in order:
   ```r
   source(here::here("scripts", "1_variability_analysis.R"))
   source(here::here("scripts", "2_spearman_correlation.R"))
   ```

Script 1 fits one Bayesian dose-response curve per site per insecticide via
MCMC. This is computationally intensive and may take a considerable time;
it uses multiple CPU cores. Re-running is cheap — curves whose fit already
exists are skipped (checkpointing). Set `skip_existing <- FALSE` at the top
of script 1 to force a full refit.

## Data format

Each CSV has one row per dose × replicate, with columns:

| Column        | Meaning                                                        |
|---------------|----------------------------------------------------------------|
| `Site`        | Anonymised laboratory code (A–P; `K_star` = K\*)               |
| `Insecticide` | Treated insecticide, or the control label for that file        |
| `Dose`        | Insecticide concentration                                      |
| `Responded`   | Number of mosquitoes dead at 24 h                              |
| `Subjects`    | Number of mosquitoes tested                                    |
| `Replicate`   | Replicate identifier                                           |
| `Mortality_perc` | Observed percentage mortality                               |

Each file contains one treated insecticide plus its matched controls.
Control labels (`Acetone`, `PY_control`, `OC_Control`) are detected
case-insensitively and are reported separately, not fitted as curves.

Two transcription errors present in the original records (two full-mortality
rows where the number responding had been keyed above the number tested)
have been corrected in these files so that responding never exceeds the
number tested.

## Outputs

Script 1 writes, per insecticide, into `outputs/<insecticide>/`:
- per-curve model fits, dose-response plots and LC10/50/90/99 posteriors
- `summary/lc50_within_site.csv`, `summary/lc50_table.csv` and plots
- `summary/control_mortality.csv`
- between- and within-site comparison files and model-check tables

Script 2 writes into `outputs/spearman/`:
- `spearman_heatmap.png`
- `spearman_correlation_matrix.csv`, `spearman_pvalues.csv`
- `spearman_uncertainty.csv` (Spearman R with 95% credible interval per pair)

## Licence

Data and code are released under [CHOOSE: e.g. CC-BY-4.0 for data, MIT for code].
Please complete this section before publishing.

## Citation

If you use these data or code, please cite the associated paper:

> [AUTHORS]. [TITLE]. [JOURNAL / preprint server], [YEAR]. [DOI]

and this archive:

> [AUTHORS]. Anonymised data and code for [TITLE]. Zenodo, [YEAR]. [DOI]
