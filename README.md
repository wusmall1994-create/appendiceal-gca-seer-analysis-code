# Appendiceal goblet cell adenocarcinoma: SEER analysis code

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22082394.svg)](https://doi.org/10.5281/zenodo.22082394)

Statistical analysis code for a registry-methodology study evaluating how nested ICD-O-3 case definitions affect reported incidence, stage distribution, and survival estimates for appendiceal goblet cell adenocarcinoma (GCA).

## Public-content boundary

This repository contains **statistical code only**. It intentionally excludes:

- SEER case-level and incidence exports;
- SEER*Stat dictionaries and session files;
- manuscript files, tables populated with study results, and publication figures;
- personal information, institutional paths, credentials, and access tokens.

SEER data cannot be redistributed from this repository. Eligible users must obtain SEER Research Data independently and comply with the SEER Research Data Agreement.

## Analysis files

- `R/00_install_packages.R`: installs missing R dependencies.
- `R/01_main_analysis.R`: cohort construction, direct age standardization with Fay–Feuer intervals, Poisson-bootstrap definition contrasts, stage and coding summaries, Kaplan–Meier estimates, Cox models, proportional-hazards checks, marginal standardization, Fine–Gray competing-risk models, age/sex interaction analyses, and manuscript-figure generation.
- `R/02_joinpoint_figure.R`: reconstructs the publication trend figure from the annual-rate output and exported final models from the NCI Joinpoint Regression Software.
- `scripts/run_analysis.ps1`: Windows wrapper for dependency installation and the main R analysis.

## Revision-phase analyses

Scripts added during peer-review revision. They analyse race/ethnicity-stratified incidence and survival, assess equivalence of histology-code survival, and generate the revision figures. Python scripts require Python 3.x with numpy, pandas, scipy, and matplotlib; the R verification script requires the `survival` package.

- `revision_analysis/qc_race_incidence_export.py`: quality-control comparison of the race-stratified SEER*Stat export (input 3 below) against the original incidence export — per-histology totals, population consistency, race distribution, and per-year totals.
- `revision_analysis/race_incidence_analysis.py`: race/origin-specific age-standardized incidence (direct standardization to the 2000 US Standard Population, Fay–Feuer intervals) for case definitions D1–D3, period rate ratios (2010–2023 vs 2000–2009), and annual rates by race/origin.
- `revision_analysis/race_incidence_figure.py`: small-multiples figure of annual race/origin-specific age-standardized incidence.
- `revision_analysis/run_phase1.py`: Table 1 characteristics by race and Hispanic ethnicity; race/ethnicity-stratified 5-year overall survival (Kaplan–Meier, log-rank, univariate Cox with Breslow ties); age-specific incidence trends and rate ratios (2020–2023 vs 2000–2004); equivalence assessment of adjusted histology-code hazard ratios against a 0.80–1.25 margin with Schoenfeld power calculations.
- `revision_analysis/make_revision_figures.py`: redraws Figure 1 (aligned axes, 2019/2021 reference lines) and builds the Figure S4 subgroup forest plot from the outputs above.
- `R/03_verify_survival_race_ethnicity.R`: independent R re-implementation of the race/ethnicity survival analysis, used to verify the Python results numerically.

Run order for the revision stage: `qc_race_incidence_export.py`, then `race_incidence_analysis.py` and `run_phase1.py`, then `race_incidence_figure.py` and `make_revision_figures.py`; `R/03_verify_survival_race_ethnicity.R` can be run any time after `raw/` is populated.

## Required private inputs

Create a local `raw/` directory and place the following independently obtained SEER*Stat exports in it:

1. `gca_incidence_age_hist_year_seer17_2000_2023.csv`
2. `gca_case_listing_seer17_2000_2023.csv`
3. `gca_incidence_age_race_nhia_hist_year_seer17_2000_2023.txt` — needed only for the revision-phase scripts; a SEER*Stat rate-session export of the same appendix cases stratified by age group, race/origin recode (NHW, NHB, NHAIAN, NHAPI, Hispanic), ICD-O-3 histology, and diagnosis year, with counts and populations. Obtain it from SEER*Stat in the same way as inputs 1–2.

The scripts expect the variable labels used in the November 2025 SEER Research Data submission. Core fields include diagnosis year, age group, sex, race/ethnicity, appendix primary site, ICD-O-3 histology, diagnostic confirmation, sequence number, summary stage, survival time, vital status, and SEER cause-specific/other-cause death classifications.

No example records are supplied because even synthetic row-level data could be mistaken for redistributable SEER content.

## Case definitions

- D1: ICD-O-3 `8243/3` (strict GCA definition).
- D2: `8243/3 + 8245/3` (historical GCC/adenocarcinoid definition).
- D3: `8243/3 + 8244/3 + 8245/3` (historical expanded-spectrum sensitivity definition).

D2 and D3 are registry-definition sensitivity analyses and must not be interpreted as modern pathologic reclassification without slide review.

## Reproduction

Tested with R 4.6.1 and these package versions:

- survival 3.8-6
- ggplot2 4.0.3
- patchwork 1.3.2
- svglite 2.2.2
- cmprsk 2.2-12

On Windows PowerShell:

```powershell
Set-Location path\to\appendiceal-gca-seer-analysis-code
.\scripts\run_analysis.ps1
```

Or run directly:

```text
Rscript R/00_install_packages.R
Rscript R/01_main_analysis.R
```

Generated files are written under `outputs/formal_next_phase/` and `outputs/revision_analysis/`, which are excluded from version control.

## Joinpoint stage

Annual age-adjusted rates produced by `R/01_main_analysis.R` can be imported into NCI Joinpoint Regression Software. The study analysis used version 6.1.0.0, annual intervals, log-linear models, supplied standard errors, uncorrelated errors, grid search, 0–4 joinpoints, 4,499 permutations, and Empirical Quantile confidence intervals for the final selected model.

After exporting the final Joinpoint model estimates to:

`outputs/joinpoint/GCA_D1_D2_D3_Joinpoint.Export.Model.Estimates.txt`

run:

```text
Rscript R/02_joinpoint_figure.R
```

Joinpoint is developed by the Surveillance Research Program, National Cancer Institute. Users should follow the official software acknowledgement and citation guidance.

## Citation

Please cite release v1.1.0 using the metadata in `CITATION.cff` or the archived versions at [doi:10.5281/zenodo.22082393](https://doi.org/10.5281/zenodo.22082393) (concept DOI, always resolves to the latest version).

## Licence

The code is released under the MIT License. This licence applies only to the code in this repository and does not alter SEER data-use restrictions or the terms governing NCI Joinpoint software.

## Reproducibility scope

The repository makes the transformations and statistical models inspectable. Exact numeric reproduction additionally requires authorized access to the same SEER database submission and matching SEER*Stat exports. No raw or processed study data are archived here.
