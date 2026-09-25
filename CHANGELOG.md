# Changes from JDG v7.26

## Sharing adaptation, 25 September 2026

### Interface and organisation

- Consolidated three input workbooks plus in-script LoQs into two sheets in one workbook: Measurements and Analytes.
- Added `config.R`, a single runner, an RStudio project and a separate Stan file.
- Added explicit group pooling/selection and analyte selection. One prior row is shared across groups.
- Added strict input checks, a validation-only mode, parameter/data snapshots, failure logs and fit checkpoints.
- Removed steroid-specific plot exclusions, menstrual-field dependencies and fixed sex labels.
- Simplified figures to population estimates, participant CVP and population-CVI sensitivity; original pooled-density graphics and their winsorisation are not reproduced.
- Kept the supplied v7.26 file untouched. This is a separate adaptation.

### Changes that may affect estimates

1. **Temporal predictor:** use explicit chronological visit numbers. The supplied raw example is in random analytical order and contains `Visita`, but no `fecha_muestra`. The former fallback could treat spreadsheet order as time. Visit gaps are now retained; quantified samples are not renumbered after censoring. Re-run the full dataset to evaluate changes in trends and exclusions.
2. **Cochran critical value:** changed the F quantile from reversed degrees of freedom to `qf(alpha/k, n-1, (k-1)*(n-1), lower.tail=FALSE)`. At k=15, n=10, alpha=0.05, the corrected critical value is about 0.173684 rather than 0.296056. This affects ANOVA screening, not the Bayesian likelihood. Sample-size balancing remains the original rounded-mean approximation. Reference: https://www.itl.nist.gov/div898/software/dataplot/refman1/auxillar/cochvari.htm
3. **Sensitivity input:** reuse the main post-trend dataset in every scenario. The default scenario reuses the main fit. In v7.26, sensitivity reloaded pre-trend data.
4. **Post-trend censoring:** recheck both subject count and censoring proportion after exclusions.
5. **Input handling:** reject duplicate participant/visit keys, nonpositive concentrations, malformed results, missing sources and conflicting cell/configured LoQs. Reject unavailable trend tests in fit mode when trend screening is requested. These cases are now explicit rather than silently coerced or carried on.

### Model and reporting

- Retained the Normal setpoint / Student-t residual model, priors, CVA/CVG relative weights, CVI hierarchy, CV transformations, bootstrap design, ANOVA CVA subtraction, HBR calculation and RCV formulas from v7.26.
- Updated Stan array syntax; compile one model and reuse it across fits.
- Save `mu` and `mu_subject` draws as well as the original parameter set so convergence summaries cover the location parameters too.
- Export diagnostics for every sensitivity fit.
- Rename the misleading Bayesian `N_final` output to `N_quantified` and distinguish it from all observations used by the model.
- Do not treat historical `cv_*` parameter names as literal CVs; many represent log-scale SDs.

No full-dataset numerical equivalence claim is made. The supplied results workbook contains only four actual example rows.

## Output and quick-start update (25 September 2026)

- Shortened `README.md` to the steps Deniz needs; moved the full instructions to `docs/TECHNICAL_NOTES.md`.
- By default, only summaries and review files are exported. `OUTPUT_DETAIL <- TRUE` restores intermediate data and draw exports; `SAVE_FITS <- TRUE` restores saved Stan fits. These settings change files written, not estimates.
