# Bayesian biological variation (one measurement per sample)

For Deniz: the example is fictional. Replace its values and priors before using results.

1. Copy input/study_template.xlsx to input/my_study.xlsx. In **Measurements** replace the example row with one row per participant and visit: subject, group, chronological sample_order and one concentration column per amino acid.
2. In **Analytes** put one row per concentration column. Use exactly matching analyte names; enter the unit, the method's **LoQ** (a LoD is not automatically a LoQ), prior CVI/CVG/CVA, their weights and sources. Enter 5 for 5%, not 0.05 or an Excel 5% cell. Replace every example value.
3. In config.R set INPUT_FILE <- "input/my_study.xlsx" and MODE <- "validate"; open bv_pipeline.Rproj in RStudio and run source("run_pipeline.R"). Inspect the reported output/input_eligibility.csv, then set MODE <- "fit" and run the same command.

The ready-to-run simulated workbook is at **input/simulated_amino_acids.xlsx** and also at examples/simulated_amino_acids.xlsx. Use the full path in config.R; a bare name is ambiguous. The input folder and its example files are included in GitHub; input/my_study.xlsx and output/ are ignored by Git.

Install once: install.packages(c("readxl", "dplyr", "lme4", "rstan", "ggplot2")). RStan needs a working C++ toolchain. Fit can be slow with RUN_SENSITIVITY and RUN_ANOVA_BOOTSTRAP enabled.

The default fit writes v7.26-shaped summary_table_formatted.csv, summary_table_numeric.csv, ANOVA and parameter tables, and the original figure families (posterior densities, subject CVI, method comparison and sensitivity). Each run gets a new output folder. Optional ANOVA presentation tables: EXPORT_ANOVA_SUMMARY <- TRUE.

With one measurement per sample, CVA and CVI depend on the analytical prior. Check diagnostics.csv and the sensitivity scenarios. The adaptation has not been numerically compared with the complete steroid study; see docs/OUTPUT_MAP.md and docs/VALIDATION.md.
