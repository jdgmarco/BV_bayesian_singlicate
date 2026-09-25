# Open bv_pipeline.Rproj, then source("run_pipeline.R").
# Command line: Rscript run_pipeline.R [validate|fit] [path/to/workbook.xlsx]
local({
  if (!file.exists("config.R")) stop("Run from the extracted project folder (or open bv_pipeline.Rproj).")
  source("config.R", local = TRUE)
  args <- if (interactive()) character() else commandArgs(trailingOnly = TRUE)
  if (length(args) > 0) MODE <- args[1]
  if (length(args) > 1) INPUT_FILE <- args[2]
  if (!MODE %in% c("validate", "fit")) stop("MODE must be validate or fit.")
  deps <- c("readxl", "dplyr")
  if (MODE == "fit") deps <- c(deps, "lme4", if (RUN_BAYESIAN) "rstan", if (MAKE_PLOTS) "ggplot2")
  missing <- deps[!vapply(deps, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Missing packages: ", paste(missing, collapse = ", "), ". See README.md.")
  suppressPackageStartupMessages(library(dplyr))
  if (MODE == "fit") suppressPackageStartupMessages(library(lme4))
  source("R/input.R", local = TRUE)
  source("R/statistics.R", local = TRUE)
  source("R/plots.R", local = TRUE)
  if (!identical(PRIOR_WEIGHT_MODE, "relative")) abort("Only relative CVA/CVG weights are supported in this adaptation.")
  if (MIN_OBS_PER_SUBJECT < 2 || MIN_SUBJECTS < 2 ||
      HIGH_CENSORING_THRESHOLD < 0 || HIGH_CENSORING_THRESHOLD >= 1)
    abort("Invalid eligibility settings in config.R.")
  if (STAN_WARMUP >= STAN_ITER || STAN_WARMUP < 0 || STAN_CHAINS < 2 || STAN_CORES < 1)
    abort("Check MCMC iteration, warmup, chain and core settings.")
  if (MODE == "fit" && RUN_SENSITIVITY && !RUN_BAYESIAN)
    abort("RUN_SENSITIVITY requires RUN_BAYESIAN = TRUE.")

  input <- read_study_inputs(INPUT_FILE)
  if (!is.null(SENSITIVITY_ANALYTES) && any(!SENSITIVITY_ANALYTES %in% input$analytes$analyte))
    abort("SENSITIVITY_ANALYTES must be included in the selected Analytes sheet.")
  dir.create(OUTPUT_ROOT, recursive = TRUE, showWarnings = FALSE)
  out <- tempfile(paste0(format(Sys.time(), "%Y%m%d_%H%M%S"), "_", MODE, "_"), tmpdir = OUTPUT_ROOT)
  dir.create(out)
  csv <- function(x, name) {
    if (!is.null(x) && nrow(x) > 0) write.csv(x, file.path(out, paste0(name, ".csv")), row.names = FALSE, na = "")
  }
  if (MODE == "fit") {
    file.copy(INPUT_FILE, file.path(out, "input_used.xlsx"))
    file.copy("config.R", file.path(out, "config_used.R"))
  }
  if (OUTPUT_DETAIL) {
    config_names <- ls(pattern = "^[A-Z][A-Z0-9_]+$")
    saveRDS(mget(config_names, envir = environment()), file.path(out, "effective_config.rds"))
    csv(data.frame(setting = config_names,
                   value = vapply(mget(config_names, envir = environment()),
                                  function(x) paste(deparse(x), collapse = " "), character(1))), "effective_config")
    csv(input$analytes, "analyte_parameters")
    csv(input$long, "parsed_measurements")
  }

  keys <- unique(input$long[c("analyte", "Sex")])
  prepared <- list(); preflight <- list(); exclusion_events <- list()
  for (i in seq_len(nrow(keys))) {
    an <- keys$analyte[i]; group <- keys$Sex[i]
    d <- input$long[input$long$analyte == an & input$long$Sex == group, , drop = FALSE]
    e <- eligible_dataset(d)
    key <- sprintf("fit_%03d", i)
    prepared[[key]] <- list(data = e$data, analyte = an, group = group, reason = e$reason)
    preflight[[key]] <- data.frame(analyte = an, group = group,
       input_rows = nrow(d), missing_results = sum(is.na(d$value)),
       numeric_results_below_loq = sum(d$numeric_below_loq),
       excluded_subjects_few_samples = length(e$excluded),
       n_subjects = e$n_subjects, n_obs = e$n_obs,
       pct_censored = e$pct_censored, status = e$reason)
    if (length(e$excluded)) exclusion_events[[length(exclusion_events) + 1L]] <-
      data.frame(analyte = an, sex = group, subject = e$excluded,
                  test = "Insufficient_samples", iter = 0L)
  }
  preflight <- bind_rows(preflight)
  csv(preflight, "input_eligibility")
  csv(bind_rows(exclusion_events), "excluded_subjects")
  print(preflight)
  message("Output directory: ", normalizePath(out))
  if (MODE == "validate") {
    writeLines(c("Input validation completed. No model was fitted.",
                 "Review input_eligibility.csv, including any ineligible groups.",
                 "Set MODE = 'fit' or use Rscript run_pipeline.R fit to estimate BV."),
               file.path(out, "RUN_STATUS.txt"))
  } else {
    if (!any(preflight$status == "eligible")) abort("No eligible analyte/group combination. See input_eligibility.csv.")
    set.seed(STAN_SEED)
    table_results_list <- list(); anova_results <- list(); anova_boots <- list()
    trend_reports <- list(); trend_summaries <- list(); fit_diagnostics <- list()
    sensitivity_rows <- list(); status_rows <- list(); posterior_draws <- list()
    if (SAVE_FITS) dir.create(file.path(out, "fits"))
    if (RUN_BAYESIAN) {
      rstan::rstan_options(auto_write = TRUE)
      options(mc.cores = min(STAN_CORES, STAN_CHAINS))
      model <- rstan::stan_model(file = "stan/bv_model.stan")
    }
    scenarios <- data.frame(
      label = c("Default", "CVI_tighter", "CVI_wider", "CVI_half", "CVI_double",
                "CVA_tighter", "CVA_wider", "CVA_half", "CVA_double"),
      prior_mult = c(1, 1, 1, .5, 2, 1, 1, 1, 1),
      weight = c(NA, .25, 1, NA, NA, NA, NA, NA, NA),
      cva_mult = c(1, 1, 1, 1, 1, 1, 1, .5, 2),
      cva_weight = c(1, 1, 1, 1, 1, .5, 2, 1, 1))
    if (OUTPUT_DETAIL) csv(scenarios, "sensitivity_scenarios")
    q <- function(x) unname(quantile(x, c(.025, .5, .975)))
    add_status <- function(an, group, stage, status, detail = "") {
      status_rows[[length(status_rows) + 1L]] <<-
        data.frame(analyte = an, group = group, stage = stage, status = status, detail = detail)
      csv(bind_rows(status_rows), "analysis_status")
    }
    record_diag <- function(fit, an, group, scenario) {
      di <- diagnose_fit(fit, paste(an, group, scenario))
      fit_diagnostics[[length(fit_diagnostics) + 1L]] <<-
        cbind(data.frame(analyte = an, group = group, scenario = scenario), as.data.frame(di))
      csv(bind_rows(fit_diagnostics), "diagnostics")
    }
    for (key in names(prepared)) {
      item <- prepared[[key]]; an <- item$analyte; group <- item$group
      if (item$reason != "eligible") {
        add_status(an, group, "eligibility", "skipped", item$reason); next
      }
      d <- item$data
      if (RUN_TREND) {
        tr <- tryCatch(trend_analysis(d, an, group), error = function(e) e)
        if (inherits(tr, "error")) {
          add_status(an, group, "trend", "failed", conditionMessage(tr)); next
        }
        trend_reports[[key]] <- tr$per_subject; trend_summaries[[key]] <- tr$summary
        if (is.na(tr$summary$global_p)) {
          add_status(an, group, "trend", "failed", "Global trend could not be assessed."); next
        }
        if (TREND_EXCLUDE && isTRUE(tr$summary$global_significant)) {
          bad <- tr$per_subject$subject[which(tr$per_subject$significant)]
          if (length(bad)) {
            d <- d[!d$subject %in% bad, , drop = FALSE]
            exclusion_events[[length(exclusion_events) + 1L]] <-
              data.frame(analyte = an, sex = group, subject = bad, test = "Trend", iter = 0L)
          }
        }
      }
      if (length(unique(d$subject)) < MIN_SUBJECTS || mean(d$censored) > HIGH_CENSORING_THRESHOLD) {
        add_status(an, group, "post_trend", "skipped", "Insufficient subjects or excessive censoring after trend exclusion.")
        next
      }
      # One authoritative analysis dataset: main Bayes AND every sensitivity fit.
      if (OUTPUT_DETAIL) csv(d, paste0(key, "_bayes_data"))
      a <- input$analytes[input$analytes$analyte == an, , drop = FALSE]
      h <- as_hypers(a)
      fit <- NULL
      if (RUN_BAYESIAN) {
        fit <- tryCatch(ejecutar_modelo_bayesiano(d, model, h, a$loq), error = function(e) e)
        if (inherits(fit, "error")) {
          add_status(an, group, "Bayesian", "failed", conditionMessage(fit)); fit <- NULL
        } else {
          record_diag(fit, an, group, "Default")
          ids <- data.frame(subject_index = as.integer(as.factor(d$subject)),
                            subject_original = d$subject, sex = d$Sex) %>% distinct() %>% arrange(subject_index)
          qs <- unname(quantile(d$value, c(.25, .5, .75)))
          table_results_list[[key]] <- list(fit = fit, sujetos_originales = ids,
            analyte = an, sex = group, n_obs = nrow(d), n_subj = length(unique(d$subject)),
            n_censored = sum(d$censored), n_quant = sum(d$censored == 0),
            conc_median = qs[2], conc_q1 = qs[1], conc_q3 = qs[3], cva_prior = a$prior_cva_pct)
          if (SAVE_FITS) saveRDS(table_results_list[[key]], file.path(out, "fits", paste0(key, "_main.rds")))
          posterior_draws[[key]] <- cbind(data.frame(analyte = an, group = group), compute_draws_table(fit))
          add_status(an, group, "Bayesian", "fitted")
        }
      }
      if (RUN_ANOVA) {
        da <- d[d$censored == 0, , drop = FALSE]
        if (RUN_OUTLIERS_ANOVA) {
          trimmed <- cochran_reed_trim(da, an, group)
          da <- trimmed$data
          if (nrow(trimmed$removed)) exclusion_events[[length(exclusion_events) + 1L]] <- trimmed$removed
        }
        if (OUTPUT_DETAIL) csv(da, paste0(key, "_anova_data"))
        ap <- calcular_anova_cv(da, cv_a_pct = a$prior_cva_pct, data_full = d)
        ap$analyte <- an; ap$group <- group; ap$n_subjects <- length(unique(da$subject))
        ap$cva_prior_pct <- a$prior_cva_pct
        anova_results[[key]] <- ap
        add_status(an, group, "ANOVA", if (any(is.finite(ap$cv_within_raw))) "fitted" else "failed",
                   if (any(is.finite(ap$cv_within_raw))) "" else "No estimable ANOVA residual variance.")
        if (RUN_ANOVA_BOOTSTRAP && any(is.finite(ap$cv_within_raw))) {
          ab <- bootstrap_anova_cv(da, cv_a_pct = a$prior_cva_pct)
          if (!is.null(ab) && nrow(ab)) {
            ab$analyte <- an; ab$group <- group; anova_boots[[key]] <- ab
          }
        }
      }
      if (RUN_SENSITIVITY && !is.null(fit) &&
          (is.null(SENSITIVITY_ANALYTES) || an %in% SENSITIVITY_ANALYTES)) {
        for (k in seq_len(nrow(scenarios))) {
          sc <- scenarios[k, ]
          sf <- if (k == 1L) fit else tryCatch(
            ejecutar_modelo_bayesiano(d, model, h, a$loq,
              weight_override = if (is.na(sc$weight)) NULL else sc$weight,
              prior_mult_override = sc$prior_mult, cva_mult_override = sc$cva_mult,
              cva_weight_mult = sc$cva_weight), error = function(e) e)
          if (inherits(sf, "error")) {
            add_status(an, group, paste0("Sensitivity:", sc$label), "failed", conditionMessage(sf)); next
          }
          if (k != 1L) record_diag(sf, an, group, sc$label)
          sq <- q(compute_draws_table(sf)$CVI_est_pct)
          sensitivity_rows[[length(sensitivity_rows) + 1L]] <- data.frame(
            analyte = an, group = group, scenario = sc$label,
            cvi_q025_pct = sq[1], cvi_median_pct = sq[2], cvi_q975_pct = sq[3],
            n_subjects = length(unique(d$subject)), n_obs = nrow(d))
          if (SAVE_FITS && k != 1L) saveRDS(sf, file.path(out, "fits", paste0(key, "_", sc$label, ".rds")))
          add_status(an, group, paste0("Sensitivity:", sc$label), "fitted")
          csv(bind_rows(sensitivity_rows), "sensitivity_population_cvi")
        }
      }
      csv(bind_rows(exclusion_events), "excluded_subjects")
    }
    csv(bind_rows(trend_summaries), "trend_summary")
    csv(bind_rows(trend_reports), "trend_by_subject")
    csv(bind_rows(exclusion_events), "excluded_subjects")
    anova <- summarise_anova(anova_results, anova_boots)
    csv(anova, "anova_summary")
    if (OUTPUT_DETAIL) csv(bind_rows(anova_boots), "anova_bootstrap_draws_fraction")
    if (length(table_results_list)) {
      formatted <- build_summary_table_formatted(table_results_list, as_hypers(input$analytes))
      numeric <- build_summary_table_numeric(table_results_list)
      for (nm in c("formatted", "numeric")) {
        tab <- get(nm)
        names(tab)[names(tab) == "Sex"] <- "Group"
        names(tab)[names(tab) == "N_final"] <- "N_quantified"
        names(tab)[names(tab) == "N_obs"] <- "N_used_in_model"
        names(tab)[names(tab) == "II_Harris_Boyd"] <- "II"
        tab$Unit <- input$analytes$unit[match(tab$Analyte, input$analytes$analyte)]
        csv(tab, paste0("bayesian_summary_", nm))
      }
      if (OUTPUT_DETAIL) csv(bind_rows(posterior_draws), "bayesian_population_draws_pct")
      subject_summary <- summarise_subjects(table_results_list)
      csv(subject_summary, "subject_cvp_pct")
      if (OUTPUT_DETAIL) {
        predictions <- bind_rows(lapply(table_results_list, function(r) {
          m <- as.data.frame(r$fit)
          bind_rows(lapply(c("dCVP_20", "dCVP_50", "dCVP_80"), function(parameter) {
            qs <- q(100 * m[[parameter]])
            data.frame(analyte = r$analyte, group = r$sex, parameter = parameter,
                       q025_pct = qs[1], median_pct = qs[2], q975_pct = qs[3])
          }))
        }))
        csv(predictions, "new_participant_cvp_percentiles")
      }
    } else subject_summary <- data.frame()
    if (MAKE_PLOTS) make_plots(out, bind_rows(posterior_draws), anova,
                              subject_summary, bind_rows(sensitivity_rows))
    status <- bind_rows(status_rows)
    failures <- sum(status$status == "failed")
    writeLines(c(paste("Fit run finished. Failed stages:", failures),
                 "Review analysis_status.csv and diagnostics.csv before interpreting estimates.",
                 "A completed fit is not a statement of adequate convergence or model validity."),
               file.path(out, "RUN_STATUS.txt"))
    if (failures) warning(failures, " stage(s) failed. See analysis_status.csv.")
  }
  if (MODE == "fit") writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
  message("Finished. Files: ", normalizePath(out))
})
