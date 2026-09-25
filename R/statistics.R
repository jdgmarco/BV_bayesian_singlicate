# Core estimators adapted from Jorge Diaz-Garzon v7.26 (August 2026).
# This file is implementation code. Configure the study in config.R and Excel.
# See CHANGELOG.md for the limited corrections and reporting changes.
# Historical cv_* Stan names often denote log-scale SDs, not linear CVs.

trend_analysis <- function(data, analyte = NA_character_, sexo = NA_character_,
                           alpha = TREND_ALPHA, min_obs = TREND_MIN_OBS,
                           global_first = TREND_GLOBAL_FIRST) {
  d <- data %>%
    filter(censored == 0, !is.na(value), value > 0) %>%
    arrange(subject, sample_order)

  empty_subj <- tibble::tibble(subject = character(), n = integer(),
                               slope_log = numeric(), p_value = numeric(),
                               significant = logical(),
                               analyte = character(), sex = character())

  # ---- 1) GLOBAL trend via lmer fixed slope -------------------------------
  global_slope <- NA_real_; global_p <- NA_real_; global_sig <- NA
  n_subj <- length(unique(d$subject))
  if (nrow(d) >= (n_subj + 3) && n_subj >= 2 &&
      length(unique(d$sample_order)) >= 3) {
    gfit <- tryCatch(
      lmer(log(value) ~ sample_order + (1 | subject), data = d, REML = FALSE),
      error = function(e) NULL)
    if (!is.null(gfit)) {
      co <- summary(gfit)$coefficients
      if ("sample_order" %in% rownames(co)) {
        global_slope <- co["sample_order", "Estimate"]
        # t value -> two-sided p via normal approx (lmerTest not assumed)
        tval <- co["sample_order", "t value"]
        global_p <- 2 * pnorm(abs(tval), lower.tail = FALSE)
        global_sig <- global_p < alpha
      }
    }
  }

  if (!is.na(global_sig)) {
    cat(sprintf("   \U0001F4C8 Global trend %s|%s: slope(log)=%.4f, p=%.3f %s\n",
                analyte, sexo, global_slope, global_p,
                ifelse(global_sig, "** SIGNIFICANT", "(no linear trend detected)")))
  }

  # ---- 2) PER-SUBJECT drill-down only if global is significant ------------
  subj_tab <- empty_subj
  if (isTRUE(global_sig) || !isTRUE(global_first)) {
    subj_tab <- d %>%
      group_by(subject) %>%
      group_modify(~{
        s <- .x; n <- nrow(s)
        if (n < min_obs || length(unique(s$sample_order)) < 3) {
          return(tibble::tibble(n = n, slope_log = NA_real_,
                                p_value = NA_real_, significant = NA))
        }
        fit <- tryCatch(lm(log(value) ~ sample_order, data = s),
                        error = function(e) NULL)
        if (is.null(fit)) {
          return(tibble::tibble(n = n, slope_log = NA_real_,
                                p_value = NA_real_, significant = NA))
        }
        co <- summary(fit)$coefficients
        if (!"sample_order" %in% rownames(co)) {
          return(tibble::tibble(n = n, slope_log = NA_real_,
                                p_value = NA_real_, significant = NA))
        }
        tibble::tibble(n = n,
                       slope_log   = co["sample_order", "Estimate"],
                       p_value     = co["sample_order", "Pr(>|t|)"],
                       significant = co["sample_order", "Pr(>|t|)"] < alpha)
      }) %>%
      ungroup() %>%
      mutate(subject = as.character(subject),
             analyte = analyte, sex = sexo)
    n_tested <- sum(!is.na(subj_tab$significant))
    n_sig    <- sum(subj_tab$significant, na.rm = TRUE)
    if (n_tested > 0)
      cat(sprintf("      \u21B3 per-subject: %d/%d subjects drift (%.0f%%)\n",
                  n_sig, n_tested, 100 * n_sig / max(n_tested, 1)))
  }

  n_tested <- sum(!is.na(subj_tab$significant))
  n_sig    <- sum(subj_tab$significant, na.rm = TRUE)
  summary_row <- data.frame(
    analyte = analyte, sex = sexo,
    global_slope_log   = global_slope,
    global_p           = global_p,
    global_significant = global_sig,
    n_subjects_tested  = n_tested,
    n_subjects_trend   = n_sig,
    pct_subjects_trend = ifelse(n_tested > 0, 100 * n_sig / n_tested, NA_real_),
    stringsAsFactors = FALSE)
  list(per_subject = subj_tab, summary = summary_row)
}

# (B) Cochran test on within-subject variances. Returns offending index or NA.
cochran_outlier <- function(var_by_subject, n_per_subject, alpha = COCHRAN_ALPHA) {
  k <- length(var_by_subject)
  if (k < 3 || any(is.na(var_by_subject)) || sum(var_by_subject) <= 0) return(NA)
  g_stat <- max(var_by_subject) / sum(var_by_subject)
  n_mean <- max(2, round(mean(n_per_subject)))
  df1 <- n_mean - 1
  # Correct numerator/denominator degrees of freedom (v7.26 had these reversed).
  Fc <- qf(alpha / k, df1, df1 * (k - 1), lower.tail = FALSE)
  c_crit <- 1 / (1 + (k - 1) / Fc)
  if (is.finite(c_crit) && g_stat > c_crit) return(which.max(var_by_subject))
  NA
}

# (B') Reed range-based test on subject means. Returns offending index or NA.
reed_outlier <- function(means) {
  k <- length(means)
  if (k < 3 || any(is.na(means))) return(NA)
  o <- order(means); sm <- means[o]
  rng <- sm[k] - sm[1]
  if (rng <= 0) return(NA)
  gap_low  <- sm[2] - sm[1]
  gap_high <- sm[k] - sm[k - 1]
  if (gap_high >= rng / 3 && gap_high >= gap_low) return(o[k])
  if (gap_low  >= rng / 3)                        return(o[1])
  NA
}

# (B'') Iterative Cochran + Reed trimming on the log scale (ANOVA only).
cochran_reed_trim <- function(d, analyte = NA, sexo = NA,
                              max_iter = OUTLIER_MAX_ITER, use_reed = REED_ENABLE) {
  removed <- list()
  work <- d %>% filter(censored == 0, !is.na(value), value > 0)
  for (iter in seq_len(max_iter)) {
    stats <- work %>% group_by(subject) %>%
      summarise(m = mean(log(value)), v = var(log(value)),
                n = dplyr::n(), .groups = "drop")
    if (nrow(stats) < 3) break
    ci <- cochran_outlier(stats$v, stats$n)
    if (!is.na(ci)) {
      bad <- stats$subject[ci]
      removed[[length(removed)+1]] <- data.frame(analyte=analyte, sex=sexo,
        iter=iter, test="Cochran", subject=as.character(bad), stringsAsFactors=FALSE)
      work <- work %>% filter(subject != bad); next
    }
    if (use_reed) {
      ri <- reed_outlier(stats$m)
      if (!is.na(ri)) {
        bad <- stats$subject[ri]
        removed[[length(removed)+1]] <- data.frame(analyte=analyte, sex=sexo,
          iter=iter, test="Reed", subject=as.character(bad), stringsAsFactors=FALSE)
        work <- work %>% filter(subject != bad); next
      }
    }
    break
  }
  removed_df <- if (length(removed)) do.call(rbind, removed) else
    data.frame(analyte=character(), sex=character(), iter=integer(),
               test=character(), subject=character())
  if (nrow(removed_df) > 0)
    cat(sprintf("   \U0001F9F9 ANOVA outliers removed: %d (Cochran/Reed)\n",
                nrow(removed_df)))
  list(data = work, removed = removed_df)
}
calcular_anova_cv <- function(data, cv_a_pct = NA_real_, data_full = data) {
  data_uncens <- data %>% filter(censored == 0)
  # [v7.17] concentration median [Q1-Q3] from FULL data (incl. imputed <LoQ),
  # computed per sex, so it matches the Bayesian table.
  conc_stats <- function(df) {
    v <- df$value[is.finite(df$value)]
    if (length(v) == 0) return(c(med = NA_real_, q1 = NA_real_, q3 = NA_real_))
    q <- stats::quantile(v, c(0.25, 0.5, 0.75), names = FALSE, na.rm = TRUE)
    c(med = q[2], q1 = q[1], q3 = q[3])
  }
  na_result <- list(
    cv_within_raw       = NA_real_,
    cv_within_corrected = NA_real_,
    cv_between          = NA_real_
  , n_final = NA_integer_, conc_median = NA_real_, conc_q1 = NA_real_, conc_q3 = NA_real_)

  fit_lmer_safe <- function(d) {
    n_obs      <- nrow(d)
    n_subjects <- length(unique(d$subject))
    if (n_obs < n_subjects + 2 || n_obs < 5 || n_subjects < 2) {
      cat(sprintf("⚠️  ANOVA skipped: n_obs=%d, n_subj=%d\n",
                  n_obs, n_subjects))
      return(na_result)
    }
    tryCatch({
      modelo <- lmer(log(value) ~ (1 | subject), data = d, REML = TRUE)
      vc <- as.data.frame(VarCorr(modelo))
      sigma_w_log <- sqrt(vc$vcov[vc$grp == "Residual"])
      sigma_b_log <- sqrt(vc$vcov[vc$grp == "subject"])

      cv_within_raw <- sqrt(exp(sigma_w_log^2) - 1)
      cv_between    <- sqrt(exp(sigma_b_log^2) - 1)

      cv_within_corrected <- NA_real_
      if (!is.na(cv_a_pct)) {
        cv_a_decimal <- cv_a_pct / 100
        diff_sq <- cv_within_raw^2 - cv_a_decimal^2
        if (diff_sq > 0) {
          cv_within_corrected <- sqrt(diff_sq)
        } else {
          cat(sprintf("  ⚠️  CV_A (%.1f%%) >= CV_within_total (%.1f%%): ",
                      100*cv_a_decimal, 100*cv_within_raw),
              "biological CV_I undetectable (set to 0)\n", sep = "")
          cv_within_corrected <- 0
        }
      }
      list(cv_within_raw = cv_within_raw,
           cv_within_corrected = cv_within_corrected,
           cv_between = cv_between,
           n_final = nrow(d))
    }, error = function(e) {
      cat("⚠️  lmer error:", e$message, "\n"); na_result
    })
  }

  if ("Sex" %in% colnames(data_uncens)) {
    sexes <- unique(data_uncens$Sex)
    if (length(sexes) == 0) {
      return(data.frame(Sex = NA_character_,
                        cv_within_raw = NA_real_,
                        cv_within_corrected = NA_real_,
                        cv_between = NA_real_))
    }
    do.call(rbind, lapply(sexes, function(s) {
      d <- subset(data_uncens, Sex == s)
      attr(d, "analyte") <- attr(data, "analyte")
      attr(d, "sexo")    <- s
      r <- fit_lmer_safe(d)
      cs <- conc_stats(subset(data_full, Sex == s)) # FULL data incl. imputed <LoQ
      data.frame(Sex = s,
                 cv_within_raw = r$cv_within_raw,
                 cv_within_corrected = r$cv_within_corrected,
                 cv_between = r$cv_between,
                 n_final = if (!is.null(r$n_final)) r$n_final else NA_integer_,
                 conc_median = unname(cs["med"]), conc_q1 = unname(cs["q1"]), conc_q3 = unname(cs["q3"]))
    }))
  } else {
    r <- fit_lmer_safe(data_uncens)
    cs <- conc_stats(data_full)
    data.frame(cv_within_raw = r$cv_within_raw,
               cv_within_corrected = r$cv_within_corrected,
               cv_between = r$cv_between,
               n_final = if (!is.null(r$n_final)) r$n_final else NA_integer_,
               conc_median = unname(cs["med"]), conc_q1 = unname(cs["q1"]), conc_q3 = unname(cs["q3"]))
  }
}

# -----------------------------------------------------------------------------
# 4b. ANOVA cluster bootstrap  [MOD-7]
# -----------------------------------------------------------------------------
bootstrap_anova_cv <- function(data, cv_a_pct = NA_real_,
                               n_boot = ANOVA_BOOTSTRAP_N,
                               seed   = ANOVA_BOOTSTRAP_SEED) {
  # Subject-level cluster bootstrap.
  # Resample subjects with replacement; refit lmer(log(value)~(1|subject))
  # on each replicate; collect back-transformed CV_within (raw and
  # CV_A-corrected). Returns long-format draws per Sex.
  set.seed(seed)
  data_uncens <- data %>% filter(censored == 0)

  do_one_sex <- function(d, sex_label) {
    subjects <- unique(d$subject)
    n_subj <- length(subjects)
    if (n_subj < 3 || nrow(d) < 8) {
      return(data.frame(
        Sex = sex_label, draw = seq_len(n_boot),
        cv_within_raw = NA_real_, cv_within_corrected = NA_real_,
        cv_between = NA_real_
      ))
    }
    out <- vector("list", n_boot)
    for (b in seq_len(n_boot)) {
      sampled_subj <- sample(subjects, size = n_subj, replace = TRUE)
      # Rename each sampled instance so duplicate subjects act as
      # independent clusters in the random-effects fit.
      blocks <- lapply(seq_along(sampled_subj), function(k) {
        rows <- d[d$subject == sampled_subj[k], , drop = FALSE]
        rows$subject <- paste0(sampled_subj[k], "_b", k)
        rows
      })
      d_boot <- do.call(rbind, blocks)

      cv_pair <- tryCatch({
        suppressMessages(suppressWarnings({
          m  <- lmer(log(value) ~ (1 | subject), data = d_boot, REML = TRUE)
          vc <- as.data.frame(VarCorr(m))
          sw <- sqrt(vc$vcov[vc$grp == "Residual"])
          sb <- sqrt(vc$vcov[vc$grp == "subject"])
          vr <- sqrt(exp(sw^2) - 1)
          vb <- sqrt(exp(sb^2) - 1)
          c(raw = if (!is.finite(vr)) NA_real_ else vr,
            between = if (!is.finite(vb)) NA_real_ else vb)
        }))
      }, error = function(e) c(raw = NA_real_, between = NA_real_))
      cv_raw     <- unname(cv_pair["raw"])
      cv_between <- unname(cv_pair["between"])

      cv_cor <- NA_real_
      if (!is.na(cv_raw) && !is.na(cv_a_pct)) {
        diff_sq <- cv_raw^2 - (cv_a_pct/100)^2
        cv_cor  <- if (diff_sq > 0) sqrt(diff_sq) else 0
      }
      out[[b]] <- data.frame(
        Sex = sex_label, draw = b,
        cv_within_raw = cv_raw, cv_within_corrected = cv_cor,
        cv_between = cv_between
      )
    }
    do.call(rbind, out)
  }

  if ("Sex" %in% colnames(data_uncens)) {
    out <- do.call(rbind, lapply(unique(data_uncens$Sex), function(s) {
      do_one_sex(data_uncens[data_uncens$Sex == s, , drop = FALSE], s)
    }))
  } else {
    out <- do_one_sex(data_uncens, NA_character_)
  }
  out <- out[is.finite(out$cv_within_raw), , drop = FALSE]
  out
}
ejecutar_modelo_bayesiano <- function(data, modelo_stan, hypers, loq,
                                     weight_override     = NULL,
                                     prior_mult_override = NULL,
                                     cva_mult_override   = NULL,
                                     cva_weight_mult     = NULL) {
  if (!"value" %in% colnames(data)) stop("Missing 'value' column")

  hypers <- as.list(hypers[1, ])
  names(hypers) <- make.names(names(hypers))
  required <- c("prior_cv_within", "prior_cv_between", "prior_cv_method",
                "weight_cv_within", "weight_cv_inter", "weight_cva")
  if (any(vapply(hypers[required], function(x) length(x) != 1 || is.na(x), logical(1))))
    stop("Empty/missing hyperparameters")

  if (!is.null(weight_override)) {
    hypers[["weight_cv_within"]] <- weight_override
  }
  if (!is.null(prior_mult_override)) {
    hypers[["prior_cv_within"]] <-
      as.numeric(hypers[["prior_cv_within"]]) * prior_mult_override
  }
  # [v7.20 SENS-C] shift the ANALYTICAL prior centre (linear CV %, pre-transform)
  if (!is.null(cva_mult_override)) {
    hypers[["prior_cv_method"]] <-
      as.numeric(hypers[["prior_cv_method"]]) * cva_mult_override
  }

  # Weak empirical prior centre for the population homeostatic setpoint.
  # Published concentration studies may represent populations with different
  # location parameters, so each analyte x sex fit is centred on its own study
  # mean. Censored observations contribute their LoQ placeholder to this centre;
  # the very broad log-scale SD of 2 makes that approximation inconsequential
  # for the dispersion parameters of interest.
  mean_value <- mean(data$value, na.rm = TRUE)
  if (mean_value <= 0) stop("Non-positive mean")

  # Røraas et al. sample-level Student-t degrees-of-freedom prior.
  prior_nu_const    <- PRIOR_NU
  prior_nu_sd_const <- PRIOR_NU_SD

  cv_to_log_sigma <- function(cv_pct) sqrt(log(1 + (cv_pct / 100)^2))

  prior_cv_between_log <- cv_to_log_sigma(as.numeric(hypers[["prior_cv_between"]]))
  prior_cv_method_log  <- cv_to_log_sigma(as.numeric(hypers[["prior_cv_method"]]))
  prior_cv_within_log  <- cv_to_log_sigma(as.numeric(hypers[["prior_cv_within"]]))

  # [v7.21 PRIOR-1] raw weight from the hypers file, then interpreted according
  # to PRIOR_WEIGHT_MODE. Under "relative" the weight is a fraction of the prior
  # sigma (0.1 = +/-10%), which is the convention the file actually uses.
  cva_weight_raw <- if ("weight_cv_method" %in% names(hypers)) {
    as.numeric(hypers[["weight_cv_method"]])
  } else if ("weight_cva" %in% names(hypers)) {
    as.numeric(hypers[["weight_cva"]])
  } else {
    0.30
  }
  prior_cv_method_sd <- if (identical(PRIOR_WEIGHT_MODE, "absolute")) {
    cva_weight_raw
  } else {
    cva_weight_raw * prior_cv_method_log
  }

  # [v7.20 SENS-C] widen/tighten the analytical prior RELATIVE to whatever it is,
  # so the override behaves identically whether prior_cv_method_sd came from the
  # hypers file (absolute) or from the 0.30 fallback (relative).
  if (!is.null(cva_weight_mult)) {
    prior_cv_method_sd <- prior_cv_method_sd * cva_weight_mult
  }

  # [v7.21 PRIOR-1] same units fix for the between-subject prior
  weight_cv_inter_raw  <- as.numeric(hypers[["weight_cv_inter"]])
  weight_cv_inter_used <- if (identical(PRIOR_WEIGHT_MODE, "absolute")) {
    weight_cv_inter_raw
  } else {
    weight_cv_inter_raw * prior_cv_between_log
  }

  data_list <- list(
    N = nrow(data), J = length(unique(data$subject)),
    subj = as.integer(as.factor(data$subject)),
    y = log(data$value),
    censored = as.integer(data$censored),
    log_loq = if (is.na(loq)) 0 else log(loq), # unused when no observation is censored
    prior_cv_between = prior_cv_between_log,
    prior_cv_method  = prior_cv_method_log,
    prior_cv_method_sd = prior_cv_method_sd,
    weight_cv_inter  = weight_cv_inter_used,   # [v7.21 PRIOR-1] units-corrected
    prior_cv_within  = prior_cv_within_log,
    weight_cv_within = as.numeric(hypers[["weight_cv_within"]]),
    mu_mean = mean_value,
    prior_nu = prior_nu_const, prior_nu_sd = prior_nu_sd_const
  )

  cat("Hyperparameters (prior CV centres in %):\n"); print(hypers)
  # [v7.21 PRIOR-1] echo the EFFECTIVE priors actually handed to Stan, so a
  # weight/units mismatch is visible in the run log instead of silent.
  cat(sprintf(
    "Prior weights [%s]: sigma_A ~ N(%.4f, %.4f) | sigma_G ~ N(%.4f, %.4f)\n",
    PRIOR_WEIGHT_MODE, prior_cv_method_log, prior_cv_method_sd,
    prior_cv_between_log, weight_cv_inter_used))
  cat(sprintf(
    "Population-mean prior: mu ~ N(log(mean observed = %.6g), 2)\n",
    mean_value))
  cat(sprintf("Censored observations: %d / %d (%.1f%%)\n",
              sum(data$censored), nrow(data),
              100 * mean(data$censored)))

  rstan::sampling(
    object = modelo_stan, data = data_list,
    iter = STAN_ITER, warmup = STAN_WARMUP, chains = STAN_CHAINS, seed = STAN_SEED,
    control = list(adapt_delta = STAN_ADAPT_DELTA,
                   max_treedepth = STAN_MAX_TREEDEPTH),
    pars = c("mu", "mu_subject", "cv_within_subject", "cv_within_original_scale",
             "cv_between", "log_mu_cv_within", "sigma_log_cv_within",
             "mu_cv_within", "sigma_cv_within",
             "sigma_method", "nu_residual",
             "dCVP_20", "dCVP_50", "dCVP_80")
  )
}

diagnose_fit <- function(fit, label = "") {
  s <- summary(fit)$summary
  max_rhat <- max(s[, "Rhat"], na.rm = TRUE)
  min_ess  <- min(s[, "n_eff"], na.rm = TRUE)
  sampler <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
  n_divergent <- sum(vapply(sampler, function(x) {
    if ("divergent__" %in% colnames(x)) sum(x[, "divergent__"]) else 0
  }, numeric(1)))
  n_max_treedepth <- sum(vapply(sampler, function(x) {
    if ("treedepth__" %in% colnames(x)) {
      sum(x[, "treedepth__"] >= STAN_MAX_TREEDEPTH)
    } else 0
  }, numeric(1)))
  cat(sprintf(
    "📊 [%s] max Rhat = %.3f | min ESS = %.0f | divergences = %d | max-treedepth hits = %d\n",
    label, max_rhat, min_ess, n_divergent, n_max_treedepth))
  if (max_rhat > 1.01)
    cat(sprintf("⚠️  Rhat > 1.01 in %s — review convergence\n", label))
  if (min_ess < 400)
    cat(sprintf("⚠️  n_eff < 400 in %s — review efficiency\n", label))
  if (n_divergent > 0)
    cat(sprintf("⚠️  %d divergent transitions in %s\n", n_divergent, label))
  if (n_max_treedepth > 0)
    cat(sprintf("⚠️  %d max-treedepth hits in %s\n", n_max_treedepth, label))
  list(max_rhat = max_rhat, min_ess = min_ess,
       n_divergent = n_divergent,
       n_max_treedepth = n_max_treedepth)
}
log_sigma_to_linear_cv <- function(s) sqrt(exp(s^2) - 1)

# [v7.15] Harris-Brown heterogeneity ratio computed PER POSTERIOR DRAW from the
# per-subject CVs (cv_within_original_scale[j]). For each draw we take the SD and
# mean across subjects of their individual CV and form HBR = 100 * SD/mean. This
# yields a full posterior for HBR (median + 95% CrI), not a median-of-medians
# approximation. Threshold for homogeneity ~22.7% (Roraas 2019, ~9.7 samples/subj).
compute_hbr_draws <- function(fit) {
  d <- as.data.frame(fit)
  cols <- grep("^cv_within_original_scale\\[", names(d), value = TRUE)
  if (length(cols) < 2) return(rep(NA_real_, nrow(d)))
  M <- as.matrix(d[, cols, drop = FALSE])           # draws x subjects
  mu <- rowMeans(M)
  sdv <- apply(M, 1, sd)
  100 * sdv / mu
}

compute_draws_table <- function(fit) {
  d <- as.data.frame(fit)
  CVA_pct      <- 100 * log_sigma_to_linear_cv(d$sigma_method)
  typical_sigma <- exp(d$log_mu_cv_within)
  CVI_est_pct  <- 100 * log_sigma_to_linear_cv(typical_sigma)
  CVI_pred_50  <- 100 * d$dCVP_50
  CVI_pred_20  <- 100 * d$dCVP_20
  CVI_pred_80  <- 100 * d$dCVP_80
  CVG_pct      <- 100 * log_sigma_to_linear_cv(d$cv_between)
  II           <- sqrt(CVA_pct^2 + CVI_est_pct^2) / CVG_pct
  # [v7.18] Asymmetric (log-normal) RCV from the log-scale SDs:
  #   sigma_total_log = sqrt(sigma_method^2 + sigma_bio^2)
  #   RCV_up = exp(+z*sqrt(2)*sigma_total_log) - 1 ; RCV_dn = exp(-...) - 1
  z_rcv          <- 1.96
  sigma_total_log <- sqrt(d$sigma_method^2 + typical_sigma^2)
  RCV_up_pct     <- 100 * (exp(+z_rcv * sqrt(2) * sigma_total_log) - 1)
  RCV_dn_pct     <- 100 * (exp(-z_rcv * sqrt(2) * sigma_total_log) - 1)
  # symmetric RCV kept for backward compatibility / internal use
  RCV_pct        <- z_rcv * sqrt(2) * sqrt(CVA_pct^2 + CVI_est_pct^2)
  data.frame(
    CVA_pct = CVA_pct, CVI_est_pct = CVI_est_pct,
    CVI_pred_20 = CVI_pred_20, CVI_pred_50 = CVI_pred_50,
    CVI_pred_80 = CVI_pred_80, CVG_pct = CVG_pct,
    II = II, RCV_pct = RCV_pct,
    RCV_up_pct = RCV_up_pct, RCV_dn_pct = RCV_dn_pct,
    # [v7.20 SENS-A] sigma_total_log exported: it is the scalar the likelihood
    # identifies under a singlicate design, and the only prior-sensitivity
    # metric invariant to the (non-linear) choice of RCV limit.
    sigma_total_log = sigma_total_log
  )
}

fmt_med_ci <- function(x, digits = 2) {
  sprintf("%.*f [%.*f, %.*f]",
          digits, median(x),
          digits, quantile(x, 0.025, names = FALSE),
          digits, quantile(x, 0.975, names = FALSE))
}

prior_posterior_contraction <- function(post_draws, prior_sd) {
  sd_post  <- sd(post_draws, na.rm = TRUE)
  if (!is.finite(prior_sd) || prior_sd == 0) return(NA_real_)
  1 - sd_post / prior_sd
}

build_summary_table_formatted <- function(results_list, data_hypers) {
  rows <- lapply(names(results_list), function(name) {
    R <- results_list[[name]]
    if (is.null(R$fit)) return(NULL)
    drws <- compute_draws_table(R$fit)
    diag <- diagnose_fit(R$fit, name)
    hyp <- data_hypers %>% filter(Magnitude == R$analyte)
    ppc <- NA_real_
    ppc_cva <- NA_real_
    if (nrow(hyp) > 0) {
      post_latent <- as.data.frame(R$fit)
      # Same latent scale as the fitted prior:
      # log_mu_cv_within ~ normal(log(prior_cv_within), weight_cv_within).
      ppc <- prior_posterior_contraction(
        post_latent$log_mu_cv_within,
        as.numeric(hyp$weight_cv_within[1]))
      # sigma_method ~ normal(prior_cv_method_log, prior_cv_method_sd) T[0,].
      cva_weight_raw <- if ("weight_cv_method" %in% names(hyp)) {
        as.numeric(hyp$weight_cv_method[1])
      } else if ("weight_cva" %in% names(hyp)) {
        as.numeric(hyp$weight_cva[1])
      } else 0.30
      cva_mu_log <- sqrt(log(1 + (as.numeric(hyp$prior_cv_method[1]) / 100)^2))
      cva_sd_log <- if (identical(PRIOR_WEIGHT_MODE, "absolute")) {
        cva_weight_raw
      } else {
        cva_weight_raw * cva_mu_log
      }
      ppc_cva <- prior_posterior_contraction(
        post_latent$sigma_method, cva_sd_log)
    }
    cva_prior <- if (!is.null(R$cva_prior)) R$cva_prior else NA_real_
    data.frame(
      Analyte = R$analyte, Sex = R$sex,
      N_subjects = R$n_subj,
      N_obs = R$n_obs,
      N_final = if (!is.null(R$n_quant)) R$n_quant else NA_integer_,
      Conc_median_IQR = if (!is.null(R$conc_median) && is.finite(R$conc_median))
        sprintf("%.3g [%.3g\u2013%.3g]", R$conc_median, R$conc_q1, R$conc_q3) else NA_character_,
      Pct_censored = sprintf("%.1f", 100 * R$n_censored / R$n_obs),
      CVA_prior_pct = if (is.na(cva_prior)) NA_character_ else sprintf("%.2f", cva_prior),
      CVI_estimated = fmt_med_ci(drws$CVI_est_pct, 2),
      CVG_pct = fmt_med_ci(drws$CVG_pct, 2),
      II_Harris_Boyd = fmt_med_ci(drws$II, 2),
      RCV_asym_pct = sprintf("%.1f / +%.1f",
                       median(drws$RCV_dn_pct), median(drws$RCV_up_pct)),
      HBR_pct = fmt_med_ci(compute_hbr_draws(R$fit), 1),
      PPC_CVI = sprintf("%.2f", ppc),
      # Contraction for the analytical latent SD; limited separation of CV_A
      # and CV_I is expected with one measurement per sample.
      PPC_CVA = sprintf("%.2f", ppc_cva),
      max_Rhat = sprintf("%.3f", diag$max_rhat),
      min_ESS = sprintf("%.0f", diag$min_ess),
      n_divergent = diag$n_divergent,
      n_max_treedepth = diag$n_max_treedepth,
      stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

build_summary_table_numeric <- function(results_list) {
  rows <- lapply(names(results_list), function(name) {
    R <- results_list[[name]]
    if (is.null(R$fit)) return(NULL)
    drws <- compute_draws_table(R$fit)
    qstats <- function(x) c(
      median = median(x),
      q025 = quantile(x, 0.025, names = FALSE),
      q975 = quantile(x, 0.975, names = FALSE))
    cvie <- qstats(drws$CVI_est_pct); cvg <- qstats(drws$CVG_pct)
    ii   <- qstats(drws$II); rcv <- qstats(drws$RCV_pct)
    rcv_up <- qstats(drws$RCV_up_pct); rcv_dn <- qstats(drws$RCV_dn_pct)
    hbr  <- qstats(compute_hbr_draws(R$fit))
    data.frame(
      Analyte = R$analyte, Sex = R$sex,
      N_subjects = R$n_subj, N_obs = R$n_obs,
      N_final = if (!is.null(R$n_quant)) R$n_quant else NA_integer_,
      Conc_median = if (!is.null(R$conc_median)) R$conc_median else NA_real_,
      Conc_Q1 = if (!is.null(R$conc_q1)) R$conc_q1 else NA_real_,
      Conc_Q3 = if (!is.null(R$conc_q3)) R$conc_q3 else NA_real_,
      Pct_censored = 100 * R$n_censored / R$n_obs,
      CVA_prior_pct = R$cva_prior,
      CVI_est_median = cvie["median"], CVI_est_q025 = cvie["q025"], CVI_est_q975 = cvie["q975"],
      CVG_median = cvg["median"], CVG_q025 = cvg["q025"], CVG_q975 = cvg["q975"],
      II_median = ii["median"], II_q025 = ii["q025"], II_q975 = ii["q975"],
      RCV_up_median = rcv_up["median"], RCV_up_q025 = rcv_up["q025"], RCV_up_q975 = rcv_up["q975"],
      RCV_dn_median = rcv_dn["median"], RCV_dn_q025 = rcv_dn["q025"], RCV_dn_q975 = rcv_dn["q975"],
      HBR_median = hbr["median"], HBR_q025 = hbr["q025"], HBR_q975 = hbr["q975"],
      stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}
