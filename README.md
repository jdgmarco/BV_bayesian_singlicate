# Biological variation pipeline — guide for Deniz

Use this R pipeline to estimate biological variation for amino acids or another analyte panel. It accepts **one analytical result per sample** and repeated samples per participant. The example values are fictitious.

## What to change

1. Copy `input/study_template.xlsx` to `input/my_study.xlsx`. In **Measurements**, replace the example row with your data: `subject` (participant ID), `group` (e.g. Female/Male), `sample_order` (chronological visit number), then one concentration column per analyte. Keep the headers in row 1. Use one row per participant and visit.
2. In **Analytes**, add exactly one row per analyte, with the same name as its column in Measurements. Enter its unit, **LoQ** (`loq`, in that unit), biological priors (`prior_cvi_pct`, `prior_cvg_pct`), method-specific analytical precision (`prior_cva_pct`), prior weights (`weight_cvi`, `weight_cvg`, `weight_cva`) and their sources. Enter CVs as numbers such as `5` for 5%, not as Excel percentages. **LoD is not automatically LoQ.** If no quantification threshold applies, leave `loq` blank.
3. In `config.R`, set `INPUT_FILE <- "input/my_study.xlsx"`. Select `MODE <- "validate"` first; set `MODE <- "fit"` when the input checks pass. Review the other options only if your design requires them.

The supplied Alanine row shows **format only**. Its result, LoQ, CVs, weights and source notes must be replaced with values appropriate to your analyte and assay. The one example observation is insufficient for estimation. By default, analysis needs at least **four valid visits per participant** and **five participants per analyte/group**; see `MIN_OBS_PER_SUBJECT` and `MIN_SUBJECTS` in `config.R`.

## Run

Open `bv_pipeline.Rproj` in RStudio. Install the required R packages once:

```r
install.packages(c("readxl", "dplyr", "lme4", "rstan", "ggplot2"))
```

RStan also requires a C++ toolchain. From the project directory, run:

```r
source("run_pipeline.R")
```

Or from a terminal: `Rscript run_pipeline.R validate input/my_study.xlsx`; after validation, `Rscript run_pipeline.R fit input/my_study.xlsx`. Each run writes a new folder under `output/`. Check `input_eligibility.csv` after validation, then fit diagnostics and `sensitivity_population_cvi.csv` after estimation.

The default output contains summary tables, eligibility, diagnostics, sensitivity and figures. For intermediate datasets and posterior draws, set `OUTPUT_DETAIL <- TRUE` in `config.R`. To save full Stan fits, set `SAVE_FITS <- TRUE`.

With one measurement per sample, the data alone do not separate analytical from biological variance: the **CVA prior from your method validation matters**. Review the sensitivity analysis and Stan diagnostics before interpreting CVI, CVG or RCV. Do not copy the steroid priors to amino acids.

## Files

- `input/study_template.xlsx`: Excel with one fictitious measurement and matching parameter row.
- `examples/simulated_amino_acids.xlsx`: larger fictitious dataset to try the pipeline.
- `config.R`: input path and analysis settings.
- `docs/TECHNICAL_NOTES.md`: column definitions, censoring rules, model details and output descriptions.
- `docs/VALIDATION.md`: checks performed and remaining validation.

Adapted from Jorge Díaz-Garzón's steroid BV pipeline v7.26 (August 2026), informed by Røraas et al., *Clinical Chemistry* 2019;65:995–1005. The adaptation has not been numerically reproduced against the complete athlete dataset; see `docs/VALIDATION.md`.
