# Regression check: participant mapping and percentage conversion, without Stan.
# Run from the project root: Rscript tests/check_subject_outputs.R
source("R/legacy_outputs.R")
draws <- data.frame(check.names = FALSE,
  "cv_within_original_scale[1]" = c(.1, .2, .3, .4, .5),
  "cv_within_original_scale[2]" = c(.2, .4, .6, .8, 1))
result <- list(fit = draws, analyte = "Alanine", sex = "Test group",
  sujetos_originales = data.frame(subject_index = c(2L, 1L),
                                  subject_original = c("P002", "P001")))
tab <- summarise_subjects(list(result))
actual <- unlist(tab[2, c("q025", "q25", "median", "q75", "q975")], use.names = FALSE)
expected <- unname(quantile(c(10, 20, 30, 40, 50), c(.025, .25, .5, .75, .975)))
stopifnot(identical(tab$subject, c("P002", "P001")),
          isTRUE(all.equal(tab$median, c(60, 30))),
          isTRUE(all.equal(actual, expected)),
          nrow(summarise_subjects(list())) == 0L)
message("Participant summaries: mapping, units and intervals passed.")
