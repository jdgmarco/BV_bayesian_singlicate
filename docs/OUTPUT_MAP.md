# Correspondence with Jorge's v7.26

This file records exactly which outputs were restored. The model in stan/bv_model.stan and the Bayesian/ANOVA summary calculations came from the supplied v7.26 script; changing the source data, screening or priors can still change estimates.

| v7.26 output | Output in this project | Correspondence |
|---|---|---|
| summary_table_formatted_* | summary_table_formatted.csv | Original builder and original column names, without an extra timestamp inside the dated output folder. |
| summary_table_numeric_* | summary_table_numeric.csv | Original builder and original column names. N_final counts quantified observations; N_obs includes censored observations used by the fit. |
| bayesian_results_by_sex_* | bayesian_results_by_sex.csv | Original parameter filtering and subject-ID mapping. |
| dcvp_predicted_by_sex_* | dcvp_predicted_by_sex.csv | Original parameter-level summaries for dCVP_20/50/80. |
| anova_results_by_sex_* | anova_results_by_sex.csv | Original ANOVA row schema: magnitude, sex, cv_a_used, raw/corrected CV and between-subject CV. |
| anova_bootstrap_summary_* | anova_bootstrap_summary.csv | Original fraction-scale CV fields and percentile summary. |
| sensitivity_table_* | sensitivity_table.csv | Original columns and original scenario labels; the default reuses the main fit. |
| anova_summary_table_* | anova_summary_table_formatted.csv and _numeric.csv | Original ANOVA table builder, opt in with EXPORT_ANOVA_SUMMARY <- TRUE as in v7.26. |
| method_comparison_* | method_comparison.csv | Bayesian and bootstrap distribution quantiles; generated with PLOT_METHOD_CMP <- TRUE. |

The six original figure files are posterior_CV_densities.png, CVi_by_subject_summary.png, method_comparison.png, sensitivity_plot.png, sensitivity_plot_reordered.png, and one CVi_per_subject_*.png per fitted analyte under plots_by_subject/. This version removes steroid-specific exclusions and captions, uses arbitrary group labels, and retains the original figure types and density/boxplot data. Its styling and display-tail handling are not pixel-identical to v7.26. Figure cropping affects only displayed draws, never tables or fits. Sensitivity figures require RUN_SENSITIVITY <- TRUE; method comparison requires a Bayesian fit and a usable ANOVA bootstrap.

Additional audit files include input_eligibility.csv, analysis_status.csv, diagnostics.csv, input_used.xlsx, config_used.R and sessionInfo.txt. Intermediate draws and participant quantile tables are exported only with OUTPUT_DETAIL <- TRUE. Trend and outlier logs are opt in through EXPORT_TREND and EXPORT_OUTLIERS. No figure is generated in MODE <- "validate".

**Differences that can change results:** the visit predictor now uses sample_order rather than spreadsheet row order; the Cochran critical value uses corrected F degrees of freedom; eligibility is rechecked after trend exclusion; the same post-trend data go into the default and all sensitivity scenarios; malformed inputs are rejected. The default sensitivity fit is reused rather than resampled, removing a pure Monte Carlo difference in that row. See CHANGELOG.md for details. The full steroid dataset was not available for numerical comparison and R/Stan was not available in the preparation environment. No claim of identical estimates or completed end-to-end figure testing is made.
