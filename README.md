# Case-definition sensitivity in cancer registry trends

Statistical analysis code for GCA nested-definition comparisons, Klatskin redistribution, simulation, and Joinpoint export audits. This repository contains source code and software documentation only. It includes no manuscript, empirical or simulated study dataset, model export, figure, contact detail, credential or local machine path. Obtain SEER inputs independently under the SEER Research Data Agreement.

## Analysis

The `methods/` directory contains the analysis scripts. Tested environment: R 4.6.1 (base R), Python 3 (standard library for the export audit), and independently installed NCI Joinpoint 6.1.0. The scripts neither install nor distribute Joinpoint. They generate simulation inputs and audit exports; they do not execute official Joinpoint or substitute R screening for permutation tests.

| Script | Analysis and required private input |
|---|---|
| reproduce_gca_bootstrap.R | GCA incidence and shared-component Poisson bootstrap; aggregate morphology-by-age annual cells |
| run_methods_extension.R | Klatskin fractional redistribution and 333-scenario simulation; cleaned seven-series aggregate input |
| audit_conservation.R | Count/rate conservation for seven Klatskin series |
| audit_official_calibration.R | Matched 2,000-series official calibration exports and simulation inputs |
| run_supplemental_controls.R | Eight control scenarios and official input generation |
| audit_additional_official.R | Official control and partial-redistribution outputs |
| run_robustness_extension.R | 18 distributional scenarios and GCA start-year inputs |
| check_robustness.R | Distributional-input checks and independent screen verification |
| audit_robustness_official.R | Matched distributional and GCA start-year official exports |
| diagnose_empirical_residuals.R | Conditional DW/ACF reconstruction for 16 empirical series |
| audit_fixed_rho_official.py | Nine positive-correlation sessions and 48 selected models |

## Execution

Use a private working directory outside this repository. `Source_Data` below means your own local inputs and exports; it is not distributed. Several audit scripts write summaries alongside those inputs, so work from copies. R arguments and filenames are explicit at the start of each script.

```text
Rscript methods/reproduce_gca_bootstrap.R Source_Data/GCA_Bootstrap_Annual_Cells.csv local_results
Rscript methods/run_methods_extension.R Source_Data/KT_all_definitions_clean_1992_2023.csv local_results
Rscript methods/audit_conservation.R Source_Data
Rscript methods/run_supplemental_controls.R Source_Data
Rscript methods/run_robustness_extension.R Source_Data
Rscript methods/check_robustness.R Source_Data
Rscript methods/audit_official_calibration.R Source_Data
Rscript methods/audit_additional_official.R Source_Data
Rscript methods/audit_robustness_official.R Source_Data
Rscript methods/diagnose_empirical_residuals.R Source_Data local_results
python methods/audit_fixed_rho_official.py Source_Data/Joinpoint_Error_Sensitivity
```

Create `local_results` first. The main simulation generates 666,000 replicates and can take time. Move generated inputs into the private source directory before dependent audits. `run_robustness_extension.R` additionally requires `GCA_joinpoint_trends_source_data.csv` there. Official-export audits require separately executed Joinpoint sessions; required export names and field labels are encoded in the scripts. Control and partial-redistribution exports belong in `Source_Data/Joinpoint_Exports`.

Seeds: main simulation 20261011; calibration 20261012; controls 20261013; distributional simulation 20261014; GCA bootstrap 20260818 (10,000 draws). Fixed-correlation empirical sessions use supplied rate SEs, log transformation, annual grid 2/2/0, 0–4 joinpoints, 4,499 permutations (seed 7160), and final-model empirical-quantile intervals (5,001 resamples). Correlations 0.1/0.2/0.3 are assumptions, not estimated residual correlation.

The fixed-rho audit expects local `Cohort_Map.csv` and `Run_Plan.csv`, the three `Joinpoint_<batch>_ErrorSensitivity.txt` inputs (cohort_id, year, rate, se), and seven official text exports plus JPO per prefix. Batches are GCA2000 (3 cohorts, 72 rows), GCA2001 (3, 69), and Klatskin1992 (10, 320); prefixes end in rho010/rho020/rho030. GCA cohort IDs are D1/D2/D3. Klatskin IDs alternate ICC/ECC at fractions 0, 25, 50, 75, 100%. Script assertions document the exact plan and map fields.

## Statistical boundaries

Historical definitions are analytic alternatives, not pathology adjudication. Fractional redistribution is deterministic and assumes independent Poisson components for variance propagation. Residual diagnostics condition on selected knots and are descriptive. R BIC selection frequencies are not nominal 5% tests. A code-only release requires independently obtained inputs and official exports to reproduce empirical findings; it is not a turnkey dataset release.

## Citation and license

Released under the MIT license. `CITATION.cff` supplies software citation metadata. The archived code is available through https://doi.org/10.5281/zenodo.22082393.
