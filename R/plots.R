# Figures corresponding to the families in v7.26, for arbitrary analytes/groups.
# CV_MAX and quantile limits crop FIGURES only; they do not change estimates.
subject_posterior_draws <- function(fit, analyte, group, scenario = NA_character_) {
  d <- as.data.frame(fit)
  cols <- grep("^cv_within_original_scale\\[", names(d), value = TRUE)
  if (!length(cols)) return(data.frame())
  data.frame(analyte = analyte, group = group, scenario = scenario,
             subject_index = rep(as.integer(gsub("[^0-9]", "", cols)), each = nrow(d)),
             cv = unlist(d[cols], use.names = FALSE), stringsAsFactors = FALSE)
}

simulate_prior_cvi <- function(prior_cv_within_pct, weight_cv_within,
                               n = N_PRIOR_SAMPLES) {
  cv_to_log_sigma <- function(cv_pct) sqrt(log(1 + (cv_pct / 100)^2))
  prior_cv_within_log <- cv_to_log_sigma(prior_cv_within_pct)
  log_mu <- rnorm(n, mean = log(prior_cv_within_log), sd = weight_cv_within)
  sigma_log <- abs(rnorm(n, mean = 0, sd = weight_cv_within))
  z <- rnorm(n, 0, 1)
  sigma_subj <- exp(log_mu + sigma_log * z)
  cv_lin <- sqrt(exp(sigma_subj^2) - 1)
  cv_lin[is.finite(cv_lin)]
}

boxplot_95 <- function(x) {
  data.frame(ymin = unname(quantile(x, .025)), lower = unname(quantile(x, .25)),
             middle = median(x), upper = unname(quantile(x, .75)),
             ymax = unname(quantile(x, .975)))
}

clip_for_display <- function(d, columns, probability, margin) {
  if (!nrow(d)) return(d)
  d <- d[is.finite(d$cv) & d$cv >= 0, , drop = FALSE]
  if (!nrow(d)) return(d)
  upper <- d %>% group_by(across(all_of(columns))) %>%
    summarise(upper = min(CV_MAX, quantile(cv, probability) * margin), .groups = "drop")
  d <- left_join(d, upper, by = columns)
  d[d$cv < d$upper, , drop = FALSE]
}

make_plots <- function(out, results, analytes, boots, sensitivity_draws,
                       scenario_levels) {
  suppressPackageStartupMessages(library(ggplot2))
  problems <- character()
  save_family <- function(name, expr) {
    tryCatch(force(expr), error = function(e) {
      problems <<- c(problems, paste0(name, ": ", conditionMessage(e)))
      warning("Figure ", name, " failed: ", conditionMessage(e), call. = FALSE)
    })
  }
  observed <- data.frame()
  if (length(results) && (PLOT_POSTERIOR_DENS || PLOT_PER_SUBJECT)) {
    observed <- dplyr::bind_rows(lapply(results, function(r) {
      dr <- subject_posterior_draws(r$fit, r$analyte, r$sex)
      if (nrow(dr)) dr$subject <- r$sujetos_originales$subject_original[
        match(dr$subject_index, r$sujetos_originales$subject_index)]
      dr
    }))
    observed <- observed[is.finite(observed$cv), , drop = FALSE]
  }

  # Original family 1: participant CVI posteriors, group colours, prior density.
  if (PLOT_POSTERIOR_DENS && nrow(observed)) save_family("posterior_CV_densities", {
    visible <- clip_for_display(observed, "analyte", PANEL_UPPER_QUANTILE,
                                PANEL_UPPER_MARGIN)
    if (nrow(visible)) {
      limits <- visible %>% group_by(analyte) %>%
        summarise(upper = max(upper), .groups = "drop")
      prior <- dplyr::bind_rows(lapply(seq_len(nrow(limits)), function(i) {
        a <- analytes[analytes$analyte == limits$analyte[i], , drop = FALSE]
        v <- simulate_prior_cvi(a$prior_cvi_pct[1], a$weight_cvi[1])
        data.frame(analyte = a$analyte[1], cv = v[v < limits$upper[i]])
      }))
      anchors <- data.frame(analyte = analytes$analyte,
                            centre = analytes$prior_cvi_pct / 100)
      p <- ggplot() +
        geom_density(data = prior, aes(x = cv), fill = "grey70", colour = "grey40",
                     alpha = .35, linewidth = .6) +
        geom_density(data = visible, aes(x = cv, fill = group, colour = group),
                     alpha = .4, linewidth = 1) +
        geom_vline(data = anchors, aes(xintercept = centre), linetype = "dashed",
                   colour = "grey30", linewidth = .6) +
        facet_wrap(~analyte, scales = "free", ncol = 3) +
        expand_limits(x = 0) +
        labs(title = "Posterior subject CVI distributions by analyte and group",
             subtitle = "Grey: prior predictive density; dashed: prior centre. Upper tails cropped for display.",
             x = "Within-subject CV (fraction)", y = "Density",
             fill = "Group", colour = "Group") +
        theme_minimal(base_size = 13) + theme(legend.position = "bottom")
      ggsave(file.path(out, "posterior_CV_densities.png"), p,
             width = 11, height = 12, dpi = 300, bg = "white")
    }
  })

  # Original family 2: participant medians, IQR and 95% CrI.
  if (PLOT_CVI_BY_SUBJECT && length(results)) save_family("CVi_by_subject_summary", {
    s <- summarise_subjects(results)
    s <- s[is.finite(s$median), , drop = FALSE]
    if (nrow(s)) {
      s <- s[order(s$analyte, s$group, s$median), , drop = FALSE]
      s$key <- factor(paste(s$analyte, s$group, s$subject, sep = " | "),
                      levels = unique(paste(s$analyte, s$group, s$subject, sep = " | ")))
      p <- ggplot(s, aes(y = key, x = median, colour = group)) +
        geom_segment(aes(x = q25, xend = q75, yend = key), linewidth = 1) +
        geom_segment(aes(x = q025, xend = q975, yend = key),
                     linewidth = .7, linetype = "dotted") +
        geom_point(size = 2.2) +
        facet_wrap(~analyte, scales = "free", ncol = 4) +
        labs(title = "Subject-specific within-subject biological variation",
             subtitle = "Posterior median CVI; IQR (solid) and 95% CrI (dotted).",
             x = "Subject CVI (%)", y = "Observed participant", colour = "Group") +
        theme_minimal(base_size = 12) +
        theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
              panel.grid.major.y = element_blank(), legend.position = "top")
      ggsave(file.path(out, "CVi_by_subject_summary.png"), p,
             width = 14, height = 10, dpi = 300, bg = "white")
    }
  })

  # Original family 3: one 95%-whisker posterior boxplot per analyte.
  if (PLOT_PER_SUBJECT && nrow(observed)) save_family("CVi_per_subject", {
    visible <- observed[observed$cv >= 0 & observed$cv < CV_MAX, , drop = FALSE]
    if (nrow(visible)) {
      plots_dir <- file.path(out, "plots_by_subject")
      dir.create(plots_dir, showWarnings = FALSE)
      for (an in unique(visible$analyte)) {
        d <- visible[visible$analyte == an, , drop = FALSE]
        d$key <- paste(d$group, d$subject, sep = " · ")
        med <- d %>% group_by(group, key) %>%
          summarise(median_cv = median(cv), .groups = "drop") %>%
          arrange(group, median_cv)
        d$key <- factor(d$key, levels = med$key)
        p <- ggplot(d, aes(x = key, y = cv, fill = group, colour = group)) +
          stat_summary(fun.data = boxplot_95, geom = "boxplot",
                       width = .6, alpha = .5, linewidth = .5) +
          facet_wrap(~group, scales = "free_x") +
          labs(title = paste("Per-subject CVI posterior —", an),
               subtitle = "Box: IQR; whiskers: 95% CrI; line: median. Upper tail cropped for display.",
               x = "Subject", y = "Within-subject CV (fraction)") +
          theme_minimal(base_size = 12) +
          theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 7),
                legend.position = "none")
        filename <- paste0("CVi_per_subject_", gsub("[^A-Za-z0-9_-]", "_", an), ".png")
        ggsave(file.path(plots_dir, filename), p,
               width = 12, height = 6, dpi = 300, bg = "white")
      }
    }
  })

  # Original family 4: three methods as distributions plus their summary CSV.
  if (PLOT_METHOD_CMP && length(results) && length(boots)) save_family("method_comparison", {
    bayes <- dplyr::bind_rows(lapply(results, function(r) {
      data.frame(Analyte = r$analyte, Sex = r$sex, method = "Bayesian",
                 cvi = compute_draws_table(r$fit)$CVI_est_pct)
    }))
    boot <- dplyr::bind_rows(boots)
    raw <- data.frame(Analyte = boot$magnitude, Sex = boot$sex_label,
                      method = "ANOVA_raw", cvi = 100 * boot$cv_within_raw)
    cor <- data.frame(Analyte = boot$magnitude, Sex = boot$sex_label,
                      method = "ANOVA_corrected", cvi = 100 * boot$cv_within_corrected)
    all <- dplyr::bind_rows(bayes, raw, cor)
    all <- all[is.finite(all$cvi), , drop = FALSE]
    if (nrow(all)) {
      all$method <- factor(all$method,
                           levels = c("Bayesian", "ANOVA_raw", "ANOVA_corrected"))
      summary <- all %>% group_by(Analyte, Sex, method) %>%
        summarise(cvi_median = median(cvi),
                  cvi_q025 = quantile(cvi, .025), cvi_q25 = quantile(cvi, .25),
                  cvi_q75 = quantile(cvi, .75), cvi_q975 = quantile(cvi, .975),
                  n_draws = n(), .groups = "drop")
      write.csv(summary, file.path(out, "method_comparison.csv"), row.names = FALSE)
      upper <- min(300, unname(quantile(all$cvi, PLOT_METHOD_CMP_QUANTILE)) *
                         PLOT_METHOD_CMP_MARGIN)
      colours <- c(Bayesian = "#1f77b4", ANOVA_raw = "#bdbdbd",
                   ANOVA_corrected = "#2ca02c")
      p <- ggplot(all, aes(x = Analyte, y = cvi, fill = method, colour = method)) +
        stat_summary(fun.data = boxplot_95, geom = "boxplot",
                     position = position_dodge(width = .75), width = .6,
                     alpha = .55, linewidth = .5) +
        facet_wrap(~Sex, ncol = 1) + coord_cartesian(ylim = c(0, upper)) +
        scale_fill_manual(values = colours) + scale_colour_manual(values = colours) +
        labs(title = "Method comparison: Bayesian vs ANOVA within-subject CV",
             subtitle = "Posterior and cluster-bootstrap distributions; median, IQR and 95% interval. Y-axis cropped for display; full values in CSV.",
             x = NULL, y = "Within-subject CV (%)", fill = "Method", colour = "Method") +
        theme_minimal(base_size = 12) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              legend.position = "bottom")
      ggsave(file.path(out, "method_comparison.png"), p,
             width = 13, height = 10, dpi = 300, bg = "white")
    }
  })

  # Original sensitivity density and the v7.26 reordered/highlighted version.
  if (PLOT_SENSITIVITY && nrow(sensitivity_draws)) save_family("sensitivity_plot", {
    d <- clip_for_display(sensitivity_draws, c("analyte", "group"),
                          SENS_PANEL_QUANTILE, SENS_PANEL_MARGIN)
    if (nrow(d)) {
      d$scenario <- factor(d$scenario, levels = scenario_levels)
      d$panel <- factor(paste(d$analyte, d$group, sep = " — "),
                        levels = unique(paste(d$analyte, d$group, sep = " — ")))
      colours <- c("Source-based default" = "#000000",
                   "Biological prior tight" = "#D55E00",
                   "Biological prior relaxed" = "#0072B2",
                   "Half prior CV" = "#CC79A7", "Double prior CV" = "#009E73",
                   "CVa tight" = "#E69F00", "CVa relaxed" = "#8C6BB1",
                   "Half prior CVa" = "#666666", "Double prior CVa" = "#56B4E9")
      base <- ggplot(d, aes(x = cv, colour = scenario, fill = scenario)) +
        geom_density(alpha = .22, linewidth = .9) +
        facet_wrap(~panel, scales = "free", ncol = 2) + expand_limits(x = 0) +
        scale_colour_manual(values = colours) + scale_fill_manual(values = colours) +
        labs(title = "Prior sensitivity — posterior subject CVI",
             subtitle = "Nine prior scenarios; upper tails cropped for display.",
             x = "Within-subject CV (fraction)", y = "Density",
             colour = "Scenario", fill = "Scenario") +
        theme_minimal(base_size = 13) + theme(legend.position = "bottom")
      height <- max(8, 2.2 * length(unique(d$analyte)))
      ggsave(file.path(out, "sensitivity_plot.png"), base,
             width = 12, height = height, dpi = 300, bg = "white")
      reference <- d[d$scenario == "Source-based default", , drop = FALSE]
      p <- ggplot(d, aes(x = cv, group = scenario)) +
        geom_density(aes(fill = scenario), colour = NA, alpha = .035) +
        geom_density(aes(colour = scenario), fill = NA, alpha = .82,
                     linewidth = .75) +
        geom_density(data = reference, aes(colour = scenario), fill = NA,
                     linewidth = 1.15, show.legend = FALSE) +
        facet_wrap(~panel, scales = "free", ncol = 2) + expand_limits(x = 0) +
        scale_colour_manual(values = colours) + scale_fill_manual(values = colours) +
        labs(title = "Prior sensitivity — posterior CVI distributions",
             subtitle = "Default in black; upper tails cropped for display.",
             x = "Within-subject CV (fraction)", y = "Density",
             colour = "Prior scenario", fill = NULL) +
        theme_minimal(base_size = 13) + theme(legend.position = "bottom") +
        guides(fill = "none")
      ggsave(file.path(out, "sensitivity_plot_reordered.png"), p,
             width = 12, height = height, dpi = 300, bg = "white")
    }
  })
  problems
}
