# Biological variation from serial measurements

This pipeline adapts Jorge Díaz-Garzón's steroid analysis (v7.26) for another panel, including amino acids, with **one analytical measurement per sample and repeated samples per participant**. The Bayesian model separates biological and analytical components using external analytical-imprecision priors. Its CVI estimates remain conditional on those priors.

**For Deniz:** edit an Excel workbook and `config.R`. You do not need to edit the Stan model or search through the statistical functions to replace steroid names.

## Start here

1. Extract the project and open `bv_pipeline.Rproj` in RStudio.
2. Install the required packages once:

   ```r
   install.packages(c("readxl", "dplyr", "lme4", "rstan", "ggplot2"))
   ```

   RStan also needs a working C++ toolchain. Follow the [official RStan installation guide](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started). The model uses modern `array[N]` Stan syntax.

3. Run `source("run_pipeline.R")`. The default is **validation only**, using the simulated example. Inspect `input_eligibility.csv` in the reported output folder.
4. To fit that example, set `MODE <- "fit"` in `config.R` and run the same command. For a first trial, setting `RUN_SENSITIVITY <- FALSE` and `RUN_ANOVA_BOOTSTRAP <- FALSE` reduces runtime. Restore them for the planned final analysis.
5. Copy `input/study_template.xlsx` to `input/my_study.xlsx`, fill it in, and change `INPUT_FILE` in `config.R`. Run validation before fitting.

From a terminal in the project folder:

```bash
Rscript run_pipeline.R validate examples/simulated_amino_acids.xlsx
Rscript run_pipeline.R fit input/my_study.xlsx
Rscript tests/check_inputs.R
```

Validation requires only `readxl` and `dplyr`. Model fitting additionally uses `lme4`, `rstan` and, for figures, `ggplot2`. Do not treat the four supplied athlete example rows as a dataset sufficient to fit BV.

## The input workbook

Keep the sheet names and header names below, or change the sheet names in `config.R`. Headers occupy **row 1**. Do not place titles above them. Blank rows are ignored, and leading/trailing spaces in headers are trimmed. Add a column and a matching parameter row for each additional analyte.

### Measurements

One row is one participant at one visit. Each analyte has its own result column.

| subject | group | sample_order | Alanine | Glycine |
|---|---|---:|---:|---:|
| P001 | Female | 1 | 310.2 | 102.5 |
| P001 | Female | 2 | 298.1 | <70 |
| P002 | Male | 1 | 325.7 | 117.3 |

These values illustrate the format only.

- `subject`: stable, pseudonymous participant ID. Store it as text in Excel if leading zeros matter. IDs must be unique across the whole study, including different groups.
- `group`: a fixed label for the participant, such as `Female`, `Male` or `All`. Other labels also work. Age is not used as a covariate.
- `sample_order`: the **chronological visit number**, a positive integer. Use collection order, never injection order, batch order or an arbitrary spreadsheet row number. Missing visits may leave gaps. This predictor expresses change per visit, not per elapsed day.
- Analyte columns: concentrations, in the unit specified on `Analytes`. Keep one analytical result per participant/visit/analyte. Replicated assays require a different data/model structure.
- Missing result: blank, `NA` or `N/A`. It is omitted for that analyte only.
- Left censoring: `<70`, `<=70`, `<LoQ` or `<LLOQ`. A numeric threshold in a cell must match that analyte's configured LoQ. Numeric positive results below LoQ are also classified as censored and counted in the input audit.
- Zeros, negative values, unrecognised text and mixed censoring limits are rejected. Decimal comma is accepted. Thousands separators are not.

The model supports one LoQ per analyte. Leave LoQ blank only if no censoring threshold applies; then all supplied results must be positive and quantified. No limit is invented. Per-sample limits, changing methods and concentration-dependent imprecision require model extensions.

### Analytes

One row per analyte. The `analyte` value must match its Measurements column exactly, including case.

| Column | What to enter |
|---|---|
| `analyte` | Result-column name, such as `Alanine`. |
| `unit` | Concentration unit for that column and its LoQ, such as `umol/L`. No automatic unit conversion is performed. |
| `loq` | Method-specific limit of quantification in that unit; positive, or blank if no threshold applies. A LoD is not automatically a LoQ. |
| `prior_cvi_pct` | External prior centre for within-subject biological CV. **Enter 15 for 15%, not 0.15.** |
| `prior_cvg_pct` | External prior centre for between-subject biological CV, also in percentage points. |
| `prior_cva_pct` | External analytical-imprecision CV appropriate to the method, concentrations and study measurement conditions. Enter 5 for 5%. |
| `weight_cvi` | Prior SD of the population location on the **log biological-SD scale**, also the half-Normal scale for between-participant heterogeneity in that SD. |
| `weight_cvg` | Relative multiplier of the prior centre of the between-subject log-scale SD. |
| `weight_cva` | Relative multiplier of the prior centre of the analytical log-scale SD. |
| `prior_source` | Source and justification for the biological priors; identify borrowed values explicitly. |
| `cva_source` | Source of analytical imprecision, such as method validation or a relevant independent precision study. |

All six numeric prior/weight values must be positive. They are ordinary numeric Excel cells, not cells entered as `5%`. A displayed Excel `5%` stores `0.05`, which is not the intended input for a 5% CV.

For CVA = 5% and `weight_cva = 0.10`, the prior centre is converted to `sigma_A0 = sqrt(log(1 + 0.05^2))`; its Normal prior SD is `0.10 * sigma_A0`, with the parameter constrained positive. The weight is **not a 95% interval**. `weight_cvi` has a different meaning and must not be interpreted as the same relative multiplier.

The current implementation uses **the same prior row for every group**. Duplicating an analyte row to specify sex-specific priors is rejected. Group-specific priors would require an explicit extension. Do not reuse the steroid prior values as amino-acid priors.

## Configuration you may change

| Setting in `config.R` | Use |
|---|---|
| `INPUT_FILE` | Your workbook path, relative to the project. |
| `MODE` | `validate` checks inputs; `fit` runs the selected estimators. |
| `ANALYTES`, `GROUPS` | `NULL` selects all, or provide explicit names. |
| `STRATIFY_BY_GROUP` | `TRUE` fits each group separately; `FALSE` pools selected participants. Pooling is a study-design decision. |
| `MIN_OBS_PER_SUBJECT`, `MIN_SUBJECTS` | Initial eligibility defaults: four valid observations per participant and five participants per analyte/group. Censored observations count as valid for this check. |
| `HIGH_CENSORING_THRESHOLD` | Skip groups with more than 80% censored observations by default. This is a configured study rule, not a universal validity threshold. |
| `RUN_TREND`, `TREND_EXCLUDE` | Check linear trends and optionally exclude participants carrying significant trends after a significant group trend. |
| `RUN_ANOVA`, `RUN_ANOVA_BOOTSTRAP` | Enable the frequentist comparison and its subject-cluster bootstrap. |
| `RUN_OUTLIERS_ANOVA` | Iterative Cochran/Reed subject removal for ANOVA only. |
| `RUN_SENSITIVITY`, `SENSITIVITY_ANALYTES` | Nine prior scenarios; by default all selected analytes are included. |
| `MAKE_PLOTS`, `SAVE_FITS` | Save figures and fitted R objects. |
| `OUTPUT_DETAIL` | Export intermediate datasets, raw draws and bootstrap draws (off by default). |
| `STAN_*`, `ANOVA_BOOTSTRAP_N` | Sampling and bootstrap settings. Default: four chains, 3000 iterations, 1500 warmup, 1000 bootstrap resamples. |

## Analysis and interpretation

The Bayesian model retains the v7.26 structure: Normal participant setpoints on the log-concentration scale, heterogeneous participant biological SDs, and Student-t observation residuals with an adaptive degrees-of-freedom parameter. The Student-t scale is adjusted so its log-scale marginal SD equals the specified total SD. Censored observations contribute a cumulative probability, not a measured concentration equal to LoQ.

With one measurement per sample, analytical and biological variances are not independently identified by those measurements. A defensible analytical prior and sensitivity analysis are therefore essential to interpretation. One common analytical CV is assumed for an analyte across subjects and concentrations. Reassess that assumption near LoQ or when precision changes across the measuring range.

Reported CVs are **normal-equivalent transformations** `100 * sqrt(exp(sigma^2) - 1)`. Under a log-Student-t observation model this is a reporting convention, not an exact finite arithmetic CV of the exponentiated Student-t distribution. Bayesian RCVs use the inherited normal-equivalent log-scale formula, not exact Student-t predictive change quantiles. ANOVA uses a log-normal model and the inherited quadratic CVA subtraction for its biological CV; that subtraction is an approximation to subtraction on the log-variance scale.

Trend screening uses quantified observations only, fits the supplied chronological visit number, and does not adjust for multiple testing. Failure to detect a linear trend does not establish physiological steady state. The same post-trend data are passed to the main Bayesian fit and all sensitivity scenarios. ANOVA subsequently discards censored results and applies its optional outlier screening, so the two estimators need not use identical observations.

Cochran and Reed remove **all observations from a flagged participant for that analyte/group**. Cochran's critical-value calculation uses the corrected F degrees of freedom and an approximate common sample size (rounded mean); unequal series lengths remain a limitation. The original initial four-observation eligibility rule is not reimposed after ANOVA discards censored results. Review the exported ANOVA datasets and counts, especially with substantial censoring.

Sensitivity scenarios: the default, tighter/wider biological priors (`weight_cvi` = 0.25/1.00), half/double biological prior centre, and half/double analytical prior centre or prior SD. These are one-at-a-time scenarios, not a factorial grid. If a source-specific weight is already below 0.25 or above 1.00, the inherited labels “tighter”/“wider” may not describe the ordering; compare the actual values or change the scenarios deliberately. Sensitivity figures show population CVI, with the same target as the sensitivity table.

Inspect Rhat, ESS, divergences and tree-depth hits for **every fitted scenario**. Recorded settings and diagnostics do not by themselves validate a model. Participant CVP, heterogeneity and the typical population CVI have different interpretations. The athlete-study HBR threshold of approximately 22.7% is not hardcoded as a general decision threshold.

## What is saved

Every run creates a new folder under `output/`.

Default output is kept to summaries, eligibility, exclusions, trend/status, diagnostics, sensitivity and figures. The intermediate exports named below are created only when `OUTPUT_DETAIL <- TRUE`; full fitted objects in `fits/` require `SAVE_FITS <- TRUE`. `input_used.xlsx`, `config_used.R` and `sessionInfo.txt` are kept for fit runs; validate runs write only the validation report and applicable exclusion rows. The original script emitted many more files; these switches change exports only, not the estimation.

| File | Meaning |
|---|---|
| `input_eligibility.csv` | Initial counts, missing data, censoring and group eligibility. |
| `parsed_measurements.csv` (detail) | Imported results plus censoring flags; the `value` field holds a LoQ placeholder for censored rows. |
| `excluded_subjects.csv` | Subjects removed by initial eligibility, trends, Cochran or Reed; absent when there are none. |
| `trend_summary.csv`, `trend_by_subject.csv` | Global trends and, when triggered, individual trends. |
| `analysis_status.csv` | Fitted, failed and skipped stages; failure reasons are retained. |
| `fit_###_bayes_data.csv`, `fit_###_anova_data.csv` (detail) | The actual observations sent to each estimator. |
| `bayesian_summary_numeric.csv`, `bayesian_summary_formatted.csv` | Population CVI, CVG, II, RCV and heterogeneity summaries. CVs and RCVs are percentages. |
| `subject_cvp_pct.csv` | Posterior CVP summaries for observed participants; these are not predictions for a new participant. |
| `new_participant_cvp_percentiles.csv` (detail) | Posterior summaries of the 20th, 50th and 80th percentiles of predicted CVP for a new participant under the hierarchy. These are population distribution percentiles, not concentration reference intervals. |
| `bayesian_population_draws_pct.csv` (detail) | Draw-level population summaries, CVs/RCVs in percent. `sigma_total_log` is a log-scale SD and `II` is dimensionless. |
| `anova_summary.csv` | Total within-subject CV, CVA-corrected biological CV, CVG, II and RCV, with bootstrap intervals where available. |
| `anova_bootstrap_draws_fraction.csv` (detail) | Underlying bootstrap CVs as fractions, not percentages. |
| `sensitivity_population_cvi.csv` | Population CVI under the nine scenarios, with numbers of observations and subjects. |
| `diagnostics.csv` | Rhat, effective sample size, divergences and maximum-tree-depth hits for main and sensitivity fits. |
| `fits/*.rds` | Checkpoints containing main results with participant mapping, and individual sensitivity Stan fits. |
| `input_used.xlsx`, `config_used.R`, `sessionInfo.txt` | Fit input, requested configuration and software environment. |
| `effective_config.*`, `analyte_parameters.csv`, `sensitivity_scenarios.csv` (detail) | Effective settings and method parameters. |

`N_used_in_model` includes quantified and censored observations. `N_quantified` counts only quantified observations. Concentration median/IQR fields retain the original descriptive convention of substituting LoQ for censored values; they are not censoring-adjusted concentration estimates.

Figures use arbitrary group labels and include all fitted analytes. There are no progesterone-specific exclusions or steroid-specific figure captions. Plot scales do not cap or winsorise posterior draws.

## Supplied examples and provenance

- `examples/simulated_amino_acids.xlsx`: 12 fictitious participants, six visits each, two groups and two analytes. All results, LoQs and priors are illustrative. **None are amino-acid reference estimates or method specifications.** Generated with NumPy `default_rng(20260925)`: log setpoints centred on 300 and 110 with SD 0.25; biological log SDs 0.18 and 0.23 multiplied by `exp(N(0,0.15))`; analytical log SD 0.05; Student-t residuals with 5 degrees of freedom, scaled to the total SD. Censoring at 5 and 70, respectively. Saved values are rounded to three decimals.
- `examples/athletes_format_example.xlsx`: the four supplied athlete example rows, paired with their group labels and the 11 supplied steroid prior/LoQ rows. It is a format example and will correctly fail sample-size eligibility. `Visita` became `sample_order`; the trailing space in `11-deoxycortisol ` was removed. The source prior values were preserved; their source notes were transcribed rather than independently revalidated.
- `input/study_template.xlsx`: input template with one fictitious Alanine result and one matching row of fictitious priors and method parameters. Replace both example rows with your measurements and analyte-specific values, then validate. The single example observation is insufficient for BV estimation.

Code origin: Jorge Díaz-Garzón, steroid BV pipeline v7.26 (August 2026), adapted from the Bayesian approach of Røraas et al., *Clinical Chemistry* 2019;65:995–1005. The accompanying athlete manuscript supplied for this adaptation is a working manuscript; no publication status or DOI is assigned here. Keep this provenance when reusing the code.

See `CHANGELOG.md` for changes that can affect results and `docs/VALIDATION.md` for checks completed and checks still requiring R/Stan. This adaptation is not a numerical reproduction of the manuscript until it has been rerun and compared using the complete athlete dataset.
