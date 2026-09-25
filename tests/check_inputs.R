# Run from the project root: Rscript tests/check_inputs.R
# No MCMC is run. Requires readxl and dplyr for workbook checks.
source("config.R")
source("R/input.R")
source("R/statistics.R")
must_fail <- function(expr) {
  caught <- tryCatch({force(expr); FALSE}, error = function(e) TRUE)
  if (!caught) stop("Expected rejection did not occur.")
}

x <- parse_results(c("1.2", "<0.5", "0.3", NA, "NA", "<LoQ", "1,5"), .5, "test")
stopifnot(identical(x$censored, c(0L, 1L, 1L, 0L, 0L, 1L, 0L)),
          identical(x$numeric_below_loq, c(FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE)),
          isTRUE(all.equal(x$value, c(1.2, .5, .5, NA, NA, .5, 1.5))))
must_fail(parse_results("<0.2", .5, "limit mismatch"))
must_fail(parse_results("0", .5, "zero"))
must_fail(parse_results("text", .5, "invalid result"))
must_fail(parse_results("<LoQ", NA_real_, "missing limit"))
stopifnot(parse_results("3", NA_real_, "no threshold")$censored == 0L)

# Known Cochran threshold, k=15, n=10, alpha=.05: Ccrit ~= .17368445.
# Largest variance 3 out of total 17 exceeds it; the former reversed-df
# implementation (critical value ~= .296056) would miss it.
stopifnot(cochran_outlier(c(3, rep(1, 14)), rep(10, 15)) == 1L,
          is.na(cochran_outlier(rep(1, 15), rep(10, 15))))

demo <- read_study_inputs("examples/simulated_amino_acids.xlsx")
stopifnot(nrow(demo$raw) == 72L, nrow(demo$analytes) == 2L,
          nrow(demo$long) == 144L, sum(demo$long$censored) == 13L)
for (an in demo$analytes$analyte) for (gr in unique(demo$long$Sex)) {
  d <- demo$long[demo$long$analyte == an & demo$long$Sex == gr, ]
  e <- eligible_dataset(d)
  stopifnot(e$reason == "eligible", e$n_subjects == 6L, e$n_obs == 36L)
}
athletes <- read_study_inputs("examples/athletes_format_example.xlsx")
stopifnot(nrow(athletes$raw) == 4L, nrow(athletes$analytes) == 11L,
          setequal(athletes$raw$sample_order, c(1, 4, 9, 11)),
          "11-deoxycortisol" %in% names(athletes$raw))
template <- read_study_inputs("input/study_template.xlsx")
stopifnot(nrow(template$raw) == 1L, nrow(template$analytes) == 1L,
          template$analytes$analyte == "Alanine",
          template$raw$subject == "EXAMPLE_001",
          eligible_dataset(template$long)$reason == "insufficient_subjects")
message("Input and Cochran checks passed. MCMC and full-data reproduction were not tested.")
