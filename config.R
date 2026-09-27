# EDIT THIS FILE AND YOUR INPUT WORKBOOK. All paths are relative to the project.
# Default: validate the simulated example without running any model.
INPUT_FILE <- "input/simulated_amino_acids.xlsx"
MEASUREMENTS_SHEET <- "Measurements"
ANALYTES_SHEET <- "Analytes"
MODE <- "validate"                    # "validate" or "fit"
OUTPUT_ROOT <- "output"

# NULL means all analytes / all group labels in the workbook.
ANALYTES <- NULL                      # e.g. c("Alanine", "Glycine")
GROUPS <- NULL                        # e.g. c("Female", "Male")
STRATIFY_BY_GROUP <- TRUE             # FALSE fits a pooled group called "All"

# Eligibility. Censored results count toward MIN_OBS_PER_SUBJECT.
MIN_OBS_PER_SUBJECT <- 4L
MIN_SUBJECTS <- 5L
HIGH_CENSORING_THRESHOLD <- 0.80       # skip if proportion is strictly greater

# Main analysis. Set MODE = "fit" after checking the validation files.
RUN_BAYESIAN <- TRUE
RUN_ANOVA <- TRUE
RUN_ANOVA_BOOTSTRAP <- TRUE
RUN_TREND <- TRUE
TREND_EXCLUDE <- TRUE                 # global AND individual trend required
RUN_OUTLIERS_ANOVA <- TRUE            # removes entire subjects for this analyte
RUN_SENSITIVITY <- TRUE
MAKE_PLOTS <- TRUE
SAVE_FITS <- FALSE                   # TRUE saves the large Stan fit objects
OUTPUT_DETAIL <- FALSE               # TRUE exports raw draws, bootstrap draws and per-fit data

# The figure families and core tables below correspond to v7.26. Each fit has
# its own output directory, so filenames need no second timestamp.
PLOT_POSTERIOR_DENS <- TRUE           # prior plus pooled subject CVI densities
PLOT_CVI_BY_SUBJECT <- TRUE           # summary of observed participant CVI
PLOT_PER_SUBJECT <- TRUE              # one posterior boxplot per analyte
PLOT_METHOD_CMP <- TRUE              # Bayesian versus ANOVA bootstrap
PLOT_SENSITIVITY <- TRUE             # density panels for the nine scenarios
EXPORT_ANOVA_SUMMARY <- FALSE        # optional v7.26 ANOVA presentation tables
EXPORT_TREND <- FALSE                # optional trend detail
EXPORT_OUTLIERS <- FALSE             # optional outlier log and counts

# Sensitivity uses exactly the same post-trend dataset as the main Bayes fit.
# NULL selects every analysed analyte. Nine one-at-a-time prior scenarios.
SENSITIVITY_ANALYTES <- NULL

# Advanced settings. Defaults follow v7.26 unless stated in CHANGELOG.md.
STAN_ITER <- 3000L
STAN_WARMUP <- 1500L
STAN_CHAINS <- 4L
STAN_SEED <- 1234L
STAN_CORES <- 4L
STAN_ADAPT_DELTA <- 0.95
STAN_MAX_TREEDEPTH <- 12L
PRIOR_NU <- 5
PRIOR_NU_SD <- 10
PRIOR_WEIGHT_MODE <- "relative"       # retained internally; do not change
TREND_GLOBAL_FIRST <- TRUE
TREND_ALPHA <- 0.05
TREND_MIN_OBS <- 4L
COCHRAN_ALPHA <- 0.05
REED_ENABLE <- TRUE
OUTLIER_MAX_ITER <- 10L
ANOVA_BOOTSTRAP_N <- 1000L
ANOVA_BOOTSTRAP_SEED <- 1234L
N_PRIOR_SAMPLES <- 4000L
CV_MAX <- 3                           # display only: upper plotting limit in fractional CV
PANEL_UPPER_QUANTILE <- 0.995
PANEL_UPPER_MARGIN <- 1.10
SENS_PANEL_QUANTILE <- 0.99
SENS_PANEL_MARGIN <- 1.05
PLOT_METHOD_CMP_QUANTILE <- 0.95
PLOT_METHOD_CMP_MARGIN <- 1.15
