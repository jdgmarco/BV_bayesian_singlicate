# Workbook ingestion and preflight checks. No statistical fitting in this file.
abort <- function(...) stop(paste0(...), call. = FALSE)

read_input_sheet <- function(path, sheet) {
  d <- as.data.frame(readxl::read_excel(path, sheet = sheet, col_types = "text",
                                      .name_repair = "minimal"),
                     stringsAsFactors = FALSE, check.names = FALSE)
  names(d) <- trimws(names(d))
  d[] <- lapply(d, function(x) { x <- trimws(x); x[x == ""] <- NA_character_; x })
  d <- d[rowSums(!is.na(d)) > 0, , drop = FALSE]
  empty_cols <- vapply(d, function(x) all(is.na(x)), logical(1))
  d <- d[, !(names(d) == "" & empty_cols), drop = FALSE]
  if (any(names(d) == "") || anyDuplicated(names(d)))
    abort(sheet, ": blank or duplicated column headers after trimming spaces.")
  d
}

require_columns <- function(d, columns, label) {
  missing <- setdiff(columns, names(d))
  if (length(missing)) abort(label, ": missing columns: ", paste(missing, collapse = ", "))
}

# Decimal comma is supported. Thousands separators and ambiguous strings are not.
strict_number <- function(x, label, allow_na = FALSE) {
  x <- trimws(as.character(x))
  x[x == ""] <- NA_character_
  pattern <- "^[+-]?([0-9]+([.,][0-9]*)?|[.,][0-9]+)([eE][+-]?[0-9]+)?$"
  bad <- !is.na(x) & !grepl(pattern, x)
  if (any(bad)) abort(label, ": invalid number(s): ", paste(head(unique(x[bad]), 5), collapse = ", "))
  y <- suppressWarnings(as.numeric(sub(",", ".", x, fixed = TRUE)))
  if (any(!is.na(y) & !is.finite(y)) || (!allow_na && anyNA(y)))
    abort(label, ": missing or non-finite numbers.")
  y
}

parse_results <- function(x, loq, label) {
  x <- trimws(as.character(x))
  missing <- is.na(x) | toupper(x) %in% c("", "NA", "N/A")
  censored <- !missing & grepl("^<", x)
  if (any(censored) && is.na(loq)) abort(label, ": censored results require a positive LoQ.")
  value <- rep(NA_real_, length(x))
  if (any(censored)) {
    tokens <- trimws(sub("^<\\s*=?\\s*", "", x[censored]))
    symbolic <- toupper(tokens) %in% c("LOQ", "LLOQ")
    explicit <- tokens[!symbolic]
    if (length(explicit)) {
      limits <- strict_number(explicit, paste0(label, " censoring limit"))
      if (any(abs(limits - loq) > 1e-8 * max(1, abs(loq))))
        abort(label, ": a cell censoring limit differs from the configured LoQ. ",
              "This model supports one limit per analyte; harmonise the data or extend the model.")
    }
    value[censored] <- loq
  }
  quantified <- !missing & !censored
  value[quantified] <- strict_number(x[quantified], label)
  if (any(value[quantified] <= 0))
    abort(label, ": zero/negative concentrations are invalid on the log scale. ",
          "Use a blank for missing data or an explicit <LoQ result for censoring.")
  below <- rep(FALSE, length(x))
  if (!is.na(loq)) below <- quantified & !is.na(value) & value < loq
  censored[below] <- TRUE
  value[below] <- loq
  data.frame(value = value, censored = as.integer(censored),
             numeric_below_loq = below, stringsAsFactors = FALSE)
}

read_study_inputs <- function(path) {
  if (!file.exists(path)) abort("Input workbook not found: ", path)
  raw <- read_input_sheet(path, MEASUREMENTS_SHEET)
  a <- read_input_sheet(path, ANALYTES_SHEET)
  require_columns(raw, c("subject", "group", "sample_order"), "Measurements")
  required <- c("analyte", "unit", "loq", "prior_cvi_pct", "prior_cvg_pct", "prior_cva_pct",
                "weight_cvi", "weight_cvg", "weight_cva", "prior_source", "cva_source")
  require_columns(a, required, "Analytes")
  if (!nrow(raw) || !nrow(a)) abort("Populate both input sheets before running the pipeline.")
  if (anyNA(a$analyte) || anyDuplicated(a$analyte)) abort("Analytes: one unique row per analyte is required.")
  if (!is.null(ANALYTES)) {
    if (length(setdiff(ANALYTES, a$analyte))) abort("ANALYTES includes a name absent from the Analytes sheet.")
    a <- a[a$analyte %in% ANALYTES, , drop = FALSE]
  }
  if (!nrow(a)) abort("No selected analytes.")
  require_columns(raw, a$analyte, "Measurements")
  numeric_cols <- c("prior_cvi_pct", "prior_cvg_pct", "prior_cva_pct", "weight_cvi", "weight_cvg", "weight_cva")
  for (col in numeric_cols) {
    a[[col]] <- strict_number(a[[col]], paste0("Analytes$", col))
    if (any(a[[col]] <= 0)) abort("Analytes$", col, ": all values must be >0.")
  }
  a$loq <- strict_number(a$loq, "Analytes$loq", allow_na = TRUE)
  if (any(!is.na(a$loq) & a$loq <= 0)) abort("LoQ must be positive or blank when no limit applies.")
  if (anyNA(a$unit) || anyNA(a$prior_source) || anyNA(a$cva_source))
    abort("Enter unit, prior_source and cva_source for every selected analyte.")
  if (anyNA(raw$subject) || anyNA(raw$group)) abort("Every measurement row needs subject and group.")
  if (!is.null(GROUPS)) {
    if (length(setdiff(GROUPS, unique(raw$group)))) abort("GROUPS contains a label absent from Measurements.")
    raw <- raw[raw$group %in% GROUPS, , drop = FALSE]
  }
  raw$sample_order <- strict_number(raw$sample_order, "Measurements$sample_order")
  if (any(raw$sample_order <= 0) || any(raw$sample_order %% 1 != 0))
    abort("sample_order must contain positive integer visit numbers.")
  if (anyDuplicated(raw[c("subject", "sample_order")]))
    abort("Duplicate subject + sample_order. The pipeline expects one measurement per sample.")
  group_count <- tapply(raw$group, raw$subject, function(x) length(unique(x)))
  if (any(group_count != 1)) abort("Each subject must belong to exactly one group.")
  raw$analysis_group <- if (STRATIFY_BY_GROUP) raw$group else "All"
  long <- lapply(seq_len(nrow(a)), function(i) {
    z <- parse_results(raw[[a$analyte[i]]], a$loq[i], a$analyte[i])
    data.frame(subject = raw$subject, Sex = raw$analysis_group,
                sample_order = raw$sample_order, analyte = a$analyte[i],
                raw_value = raw[[a$analyte[i]]], z, stringsAsFactors = FALSE)
  })
  list(raw = raw, analytes = a, long = dplyr::bind_rows(long))
}

as_hypers <- function(a) {
  data.frame(Magnitude = a$analyte, prior_cv_within = a$prior_cvi_pct,
             prior_cv_between = a$prior_cvg_pct, prior_cv_method = a$prior_cva_pct,
             weight_cv_within = a$weight_cvi, weight_cv_inter = a$weight_cvg,
             weight_cva = a$weight_cva, stringsAsFactors = FALSE)
}

eligible_dataset <- function(d) {
  valid <- d[!is.na(d$value), , drop = FALSE]
  counts <- table(valid$subject)
  excluded <- names(counts)[counts < MIN_OBS_PER_SUBJECT]
  absent <- setdiff(unique(d$subject), names(counts))
  excluded <- unique(c(excluded, absent))
  kept <- valid[!valid$subject %in% excluded, , drop = FALSE]
  nsub <- length(unique(kept$subject))
  cens <- if (nrow(kept)) mean(kept$censored) else NA_real_
  reason <- if (nsub < MIN_SUBJECTS) "insufficient_subjects" else
    if (cens > HIGH_CENSORING_THRESHOLD) "high_censoring" else "eligible"
  list(data = kept, excluded = excluded, reason = reason,
       n_subjects = nsub, n_obs = nrow(kept), pct_censored = 100 * cens)
}
