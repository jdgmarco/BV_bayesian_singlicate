# Validation record

Prepared 25 September 2026.

## Completed

- Read all three supplied input workbooks, the v7.26 R script and the accompanying manuscript.
- Counted four nonempty result rows, 30 participant records and 11 analyte-prior records in the supplied files. Blank formatted rows are not treated as observations.
- Matched the four result rows to participant groups and retained their chronological `Visita` values (1, 4, 9 and 11).
- Compared the prepared source and simulation matrices with the supplied workbooks. The four real athlete observations inspected during preparation have been removed from the GitHub package.
- Checked prepared workbook headers, subject/visit uniqueness, numeric prior fields and worksheet previews; the four-row athlete example is no longer distributed.
- Verified the simulated example contains 72 observations (12 participants, six visits), two analytes, six participants per group and 13 censored cells.
- Confirmed the extracted Stan model text is identical to the original statistical model apart from the modern array declaration syntax and removal of its R string wrapper.
- Independently recomputed the corrected Cochran threshold using SciPy: 0.17368445087533926 for k=15, n=10, alpha=0.05.
- Checked matching string/bracket delimiters in all R and Stan files and checked project file references. This is a static check, not an R syntax-parser or execution test.
- Inspected main/sensitivity data flow: every sensitivity scenario receives the same post-trend dataset as the main Bayesian fit, and the default scenario reuses that fit.
- Compared complete function bodies to v7.26: `build_summary_table_formatted`, `build_summary_table_numeric`, `compute_draws_table` and `calcular_anova_cv` are text-identical. The original `build_anova_summary_table` is included separately in `R/legacy_outputs.R`.
- Restored the v7.26 table names/columns and all figure families in source code. Checked that the two tracked simulated workbooks have identical bytes; these new exports/plots still require execution in R.

## Pending in an R environment

Review on 27 September 2026 found and repaired a missing participant-summary
helper, which would have stopped the fit during output generation. The regression
script `tests/check_subject_outputs.R` checks subject mapping, percentage units
and quantiles with fixed draws. It has not been executed here.

R and RStan are not installed in the preparation environment. The following have **not** been executed:

1. R syntax parsing and `Rscript tests/check_inputs.R`.
2. `Rscript run_pipeline.R validate input/simulated_amino_acids.xlsx`.
3. Stan compilation, MCMC sampling, ANOVA and bootstrap execution, and production of every restored table/figure from actual fitted objects.
4. Examination of convergence, posterior predictive adequacy and sensitivity results for a complete study.
5. Reproduction/comparison with the original athlete analysis using the full raw dataset. The supplied four-row example cannot establish numerical equivalence.

The included R checks exercise censoring, missing/invalid values, mismatched LoQs, the Cochran correction, demo workbook ingestion, and ingestion of the one-row template with its expected insufficient-subjects status. Passing them would not alone validate scientific inference. Keep the validation record updated after running the complete analysis locally.
