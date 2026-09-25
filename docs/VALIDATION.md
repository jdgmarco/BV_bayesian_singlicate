# Validation record

Prepared 25 September 2026.

## Completed

- Read all three supplied input workbooks, the v7.26 R script and the accompanying manuscript.
- Counted four nonempty result rows, 30 participant records and 11 analyte-prior records in the supplied files. Blank formatted rows are not treated as observations.
- Matched the four result rows to participant groups and retained their chronological `Visita` values (1, 4, 9 and 11).
- Compared every populated cell in each exported workbook with the prepared source/simulation matrices. Text matched exactly; numeric cells matched at relative tolerance 1e-14 and absolute tolerance 1e-15. No analytical values were rounded in storage for the athlete example.
- Checked exported workbook headers, subject/visit uniqueness, numeric prior fields, 11 analyte-name matches, frozen headers, and all six worksheet previews.
- Verified the simulated example contains 72 observations (12 participants, six visits), two analytes, six participants per group and 13 censored cells.
- Confirmed the extracted Stan model text is identical to the original statistical model apart from the modern array declaration syntax and removal of its R string wrapper.
- Independently recomputed the corrected Cochran threshold using SciPy: 0.17368445087533926 for k=15, n=10, alpha=0.05.
- Checked matching string/bracket delimiters in all R and Stan files and checked project file references. This is a static check, not an R syntax-parser or execution test.
- Inspected main/sensitivity data flow: every sensitivity scenario receives the same post-trend dataset as the main Bayesian fit, and the default scenario reuses that fit.

## Pending in an R environment

R and RStan are not installed in the preparation environment. The following have **not** been executed:

1. R syntax parsing and `Rscript tests/check_inputs.R`.
2. `Rscript run_pipeline.R validate examples/simulated_amino_acids.xlsx`.
3. Stan compilation, MCMC sampling, ANOVA and bootstrap execution, and plotting from actual fitted objects.
4. Examination of convergence, posterior predictive adequacy and sensitivity results for a complete study.
5. Reproduction/comparison with the original athlete analysis using the full raw dataset. The supplied four-row example cannot establish numerical equivalence.

The included R checks exercise censoring, missing/invalid values, mismatched LoQs, the Cochran correction, demo workbook ingestion, and ingestion of the one-row template with its expected insufficient-subjects status. Passing them would not alone validate scientific inference. Keep the validation record updated after running the complete analysis locally.
