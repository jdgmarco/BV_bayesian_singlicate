# Export helpers retained from the v7.26 steroid pipeline.
# The ANOVA table builder below is copied verbatim from v7.26.
# The runner supplies its original field names (magnitude, sex, cv_a_used, etc.).

build_anova_summary_table <- function(anova_results_list,
                                      anova_bootstrap_list = NULL,
                                      order_vec = NULL) {
  if (length(anova_results_list) == 0) return(NULL)
  pts <- dplyr::bind_rows(anova_results_list)
  # normalise analyte/sex columns
  if (!"Analyte" %in% names(pts) && "magnitude" %in% names(pts))
    pts$Analyte <- pts$magnitude
  if (!"Sex" %in% names(pts) && "sex" %in% names(pts))
    pts$Sex <- pts$sex

  # optional bootstrap CIs (raw, corrected, between)
  boot <- NULL
  if (!is.null(anova_bootstrap_list) && length(anova_bootstrap_list) > 0) {
    boot_raw <- dplyr::bind_rows(anova_bootstrap_list) %>%
      mutate(Sex = ifelse(is.na(Sex), sex_label, Sex))
    if (!"cv_between" %in% names(boot_raw)) boot_raw$cv_between <- NA_real_
    boot <- boot_raw %>%
      group_by(magnitude, Sex) %>%
      summarise(
        cvi_raw_lo = 100*quantile(cv_within_raw,       0.025, names = FALSE, na.rm = TRUE),
        cvi_raw_hi = 100*quantile(cv_within_raw,       0.975, names = FALSE, na.rm = TRUE),
        cvi_cor_lo = 100*quantile(cv_within_corrected, 0.025, names = FALSE, na.rm = TRUE),
        cvi_cor_hi = 100*quantile(cv_within_corrected, 0.975, names = FALSE, na.rm = TRUE),
        cvg_lo     = 100*quantile(cv_between, 0.025, names = FALSE, na.rm = TRUE),
        cvg_hi     = 100*quantile(cv_between, 0.975, names = FALSE, na.rm = TRUE),
        .groups = "drop") %>%
      rename(Analyte = magnitude)
  }

  num <- pts %>%
    transmute(
      Analyte, Sex,
      N_final      = if ("n_final" %in% names(.)) n_final else NA_integer_,
      Conc_median  = if ("conc_median" %in% names(.)) conc_median else NA_real_,
      Conc_Q1      = if ("conc_q1" %in% names(.)) conc_q1 else NA_real_,
      Conc_Q3      = if ("conc_q3" %in% names(.)) conc_q3 else NA_real_,
      CVA_prior_pct = cv_a_used,
      CVI_raw_pct  = 100 * cv_within_raw,
      CVI_corr_pct = 100 * cv_within_corrected,
      CVG_pct      = 100 * cv_between
    ) %>%
    mutate(
      II      = sqrt(CVA_prior_pct^2 + CVI_corr_pct^2) / CVG_pct,
      RCV_pct = sqrt(2) * 1.96 * CVI_raw_pct,
      # [v7.19] For the singlicate ANOVA, cv_within_raw already IS the total
      # intra-individual variation (CV_I + CV_A confounded; they cannot be
      # separated without replicates). Use it directly for the operative RCV,
      # rather than re-combining CV_A(prior) with CV_I(corrected), which would
      # subtract and re-add CV_A and reintroduce the prior's error twice.
      sigma_tot_log = sqrt(log(1 + (CVI_raw_pct/100)^2)),
      RCV_up_pct = 100 * (exp(+1.96 * sqrt(2) * sigma_tot_log) - 1),
      RCV_dn_pct = 100 * (exp(-1.96 * sqrt(2) * sigma_tot_log) - 1)
    )

  if (!is.null(boot)) {
    num <- num %>% left_join(boot, by = c("Analyte", "Sex"))
  } else {
    num <- num %>% mutate(cvi_raw_lo=NA_real_, cvi_raw_hi=NA_real_,
                          cvi_cor_lo=NA_real_, cvi_cor_hi=NA_real_,
                          cvg_lo=NA_real_, cvg_hi=NA_real_)
  }

  # ordering
  if (is.null(order_vec)) order_vec <- unique(num$Analyte)
  num <- num %>%
    mutate(analyte_rank = match(Analyte, order_vec),
           sex_rank = ifelse(Sex == "Female", 1, 2)) %>%
    arrange(analyte_rank, sex_rank)

  fmt   <- function(x, d = 1) ifelse(is.na(x), "\u2014", sprintf(paste0("%.", d, "f"), x))
  fmtci <- function(med, lo, hi, d = 1)
    ifelse(is.na(med), "\u2014",
           ifelse(is.na(lo) | is.na(hi), fmt(med, d),
                  sprintf("%s [%s, %s]", fmt(med, d), fmt(lo, d), fmt(hi, d))))

  formatted <- num %>%
    transmute(
      Analyte, Sex,
      N_final      = ifelse(is.na(N_final), "—", as.character(N_final)),
      Conc_median_IQR = ifelse(is.na(Conc_median), "—",
                          sprintf("%.3g [%.3g–%.3g]", Conc_median, Conc_Q1, Conc_Q3)),
      CVA_prior_pct = fmt(CVA_prior_pct, 2),
      CVI_raw_pct  = fmtci(CVI_raw_pct,  cvi_raw_lo, cvi_raw_hi),
      CVI_corr_pct = fmtci(CVI_corr_pct, cvi_cor_lo, cvi_cor_hi),
      CVG_pct      = fmtci(CVG_pct,      cvg_lo,     cvg_hi),
      II           = fmt(II, 2),
      RCV_asym_pct = ifelse(is.na(RCV_up_pct), "—",
                       sprintf("%.1f / +%.1f", RCV_dn_pct, RCV_up_pct))
    )

  list(formatted = formatted,
       numeric   = num %>% select(Analyte, Sex, N_final, Conc_median, Conc_Q1, Conc_Q3,
                                  CVA_prior_pct,
                                  CVI_raw_pct, CVI_corr_pct, CVG_pct, II,
                                  RCV_up_pct, RCV_dn_pct,
                                  cvi_raw_lo, cvi_raw_hi, cvi_cor_lo, cvi_cor_hi,
                                  cvg_lo, cvg_hi))
}

summarise_subjects <- function(results) {
  rows <- lapply(results, function(r) {
    draws <- as.data.frame(r$fit)
    ids <- r$sujetos_originales
    dplyr::bind_rows(lapply(seq_len(nrow(ids)), function(j) {
      parameter <- paste0("cv_within_original_scale[", ids$subject_index[j], "]")
      if (!parameter %in% names(draws))
        stop("Missing posterior parameter: ", parameter, call. = FALSE)
      qs <- unname(quantile(100 * draws[[parameter]], c(.025, .25, .5, .75, .975)))
      data.frame(analyte = r$analyte, group = r$sex, subject = ids$subject_original[j],
                 q025 = qs[1], q25 = qs[2], median = qs[3], q75 = qs[4], q975 = qs[5])
    }))
  })
  dplyr::bind_rows(rows)
}

# The per-parameter and new-participant Stan summaries use the same filtering,
# participant mapping and column order as v7.26.
collect_bayesian_parameter_tables <- function(results_list) {
  rows <- lapply(results_list, function(R) {
    if (is.null(R$fit)) return(NULL)
    df_all <- as.data.frame(summary(R$fit)$summary)
    df_all$parameter <- rownames(df_all)
    df_all$analyte <- R$analyte
    df_all$sex <- R$sex
    pattern_main <- paste0("cv_within_subject|cv_within_original_scale|cv_between|",
                           "log_mu_cv_within|sigma_log_cv_within|",
                           "mu_cv_within|sigma_cv_within|sigma_method|nu_residual")
    df_main <- df_all[grep(pattern_main, df_all$parameter), , drop = FALSE]
    if (nrow(df_main)) {
      df_main$subject_index <- ifelse(
        grepl("^cv_within_(subject|original_scale)\\[", df_main$parameter),
        as.numeric(gsub("[^0-9]", "", df_main$parameter)), NA)
      df_main <- merge(df_main, R$sujetos_originales, by = "subject_index",
                       all.x = TRUE, suffixes = c("", ".y"))
      if ("sex.y" %in% names(df_main)) df_main$sex.y <- NULL
      df_main <- df_main[, c("analyte", "sex", "subject_original", "parameter",
                             setdiff(names(df_main), c("analyte", "sex",
                                                       "subject_original", "parameter",
                                                       "subject_index"))), drop = FALSE]
    }
    df_dcvp <- df_all[grep("^dCVP_", df_all$parameter), , drop = FALSE]
    if (nrow(df_dcvp)) {
      df_dcvp$subject_original <- NA
      df_dcvp <- df_dcvp[, c("analyte", "sex", "subject_original", "parameter",
                             setdiff(names(df_dcvp), c("analyte", "sex",
                                                       "subject_original", "parameter"))), drop = FALSE]
    }
    list(main = df_main, predicted = df_dcvp)
  })
  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) return(list(main = data.frame(), predicted = data.frame()))
  list(main = dplyr::bind_rows(lapply(rows, `[[`, "main")),
       predicted = dplyr::bind_rows(lapply(rows, `[[`, "predicted")))
}

summarise_anova_bootstrap_v726 <- function(bootstrap_list) {
  if (!length(bootstrap_list)) return(data.frame())
  all <- dplyr::bind_rows(bootstrap_list)
  if (!"cv_between" %in% names(all)) all$cv_between <- NA_real_
  all %>%
    group_by(magnitude, Sex, sex_label) %>%
    summarise(
      n_draws = n(),
      cv_within_raw_median = median(cv_within_raw),
      cv_within_raw_q025 = quantile(cv_within_raw, .025, names = FALSE),
      cv_within_raw_q975 = quantile(cv_within_raw, .975, names = FALSE),
      cv_within_corr_median = median(cv_within_corrected, na.rm = TRUE),
      cv_within_corr_q025 = quantile(cv_within_corrected, .025, names = FALSE, na.rm = TRUE),
      cv_within_corr_q975 = quantile(cv_within_corrected, .975, names = FALSE, na.rm = TRUE),
      cv_between_median = median(cv_between, na.rm = TRUE),
      cv_between_q025 = quantile(cv_between, .025, names = FALSE, na.rm = TRUE),
      cv_between_q975 = quantile(cv_between, .975, names = FALSE, na.rm = TRUE),
      .groups = "drop"
    )
}
