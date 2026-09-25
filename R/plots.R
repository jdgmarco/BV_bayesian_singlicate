# Summary and presentation helpers. Estimates are never clipped for plotting.
summarise_anova <- function(results, boots) {
  if (!length(results)) return(data.frame())
  a <- dplyr::bind_rows(results)
  a$total_within_pct <- 100 * a$cv_within_raw
  a$cvi_corrected_pct <- 100 * a$cv_within_corrected
  a$cvg_pct <- 100 * a$cv_between
  a$II <- sqrt(a$cva_prior_pct^2 + a$cvi_corrected_pct^2) / a$cvg_pct
  s <- sqrt(log1p(a$cv_within_raw^2))
  a$rcv_down_pct <- 100 * expm1(-1.96 * sqrt(2) * s)
  a$rcv_up_pct <- 100 * expm1(1.96 * sqrt(2) * s)
  if (length(boots)) {
    b <- dplyr::bind_rows(boots) %>% dplyr::group_by(analyte, group) %>%
      dplyr::summarise(n_boot_valid = dplyr::n(),
        cvi_corrected_q025_pct = 100 * quantile(cv_within_corrected, .025, na.rm = TRUE),
        cvi_corrected_q975_pct = 100 * quantile(cv_within_corrected, .975, na.rm = TRUE),
        total_within_q025_pct = 100 * quantile(cv_within_raw, .025, na.rm = TRUE),
        total_within_q975_pct = 100 * quantile(cv_within_raw, .975, na.rm = TRUE),
        cvg_q025_pct = 100 * quantile(cv_between, .025, na.rm = TRUE),
        cvg_q975_pct = 100 * quantile(cv_between, .975, na.rm = TRUE), .groups = "drop")
    a <- dplyr::left_join(a, b, by = c("analyte", "group"))
  }
  a[, setdiff(names(a), c("Sex", "cv_within_raw", "cv_within_corrected", "cv_between")), drop = FALSE]
}

summarise_subjects <- function(results) {
  dplyr::bind_rows(lapply(results, function(r) {
    m <- as.data.frame(r$fit)
    ids <- r$sujetos_originales
    do.call(rbind, lapply(seq_len(nrow(ids)), function(j) {
      v <- 100 * m[[paste0("cv_within_original_scale[", ids$subject_index[j], "]")]]
      qs <- unname(quantile(v, c(.025, .25, .5, .75, .975)))
      data.frame(analyte = r$analyte, group = r$sex, subject = ids$subject_original[j],
                  q025 = qs[1], q25 = qs[2], median = qs[3], q75 = qs[4], q975 = qs[5])
    }))
  }))
}

make_plots <- function(out, posterior, anova, subjects, sensitivity) {
  suppressPackageStartupMessages(library(ggplot2))
  theme_set(theme_minimal(base_size = 11))
  if (nrow(posterior)) {
    b <- posterior %>% group_by(analyte, group) %>% summarise(
      median = median(CVI_est_pct), lo = quantile(CVI_est_pct, .025),
      hi = quantile(CVI_est_pct, .975), .groups = "drop")
    b$method <- "Bayesian CVI (95% CrI)"
    cmp <- b
    if (nrow(anova) && "cvi_corrected_q025_pct" %in% names(anova)) {
      a <- anova %>% transmute(analyte, group, median = cvi_corrected_pct,
                              lo = cvi_corrected_q025_pct, hi = cvi_corrected_q975_pct,
                              method = "CVA-corrected ANOVA (95% CI)")
      raw <- anova %>% transmute(analyte, group, median = total_within_pct,
                                lo = total_within_q025_pct, hi = total_within_q975_pct,
                                method = "ANOVA total within (95% CI)")
      cmp <- bind_rows(cmp, a, raw)
    }
    p <- ggplot(cmp, aes(analyte, median, colour = method)) +
      geom_pointrange(aes(ymin = lo, ymax = hi), position = position_dodge(width = .55)) +
      facet_wrap(~group, scales = "free_y") +
      labs(x = NULL, y = "CV (%)", colour = NULL, title = "Population variation estimates") +
      theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "bottom")
    ggsave(file.path(out, "population_comparison.png"), p, width = 12, height = 7, dpi = 200)
  }
  if (nrow(subjects)) {
    p <- ggplot(subjects, aes(x = median, y = factor(subject), colour = group)) +
      geom_segment(aes(x = q025, xend = q975, yend = factor(subject)), linewidth = .4) +
      geom_segment(aes(x = q25, xend = q75, yend = factor(subject)), linewidth = 1.3) +
      geom_point() + facet_wrap(~analyte + group, scales = "free", ncol = 2) +
      labs(x = "Participant CVP (%)", y = "Subject", title = "Observed participants",
           subtitle = "Posterior median, IQR (thick) and 95% CrI (thin)") +
      theme(legend.position = "none")
    ggsave(file.path(out, "subject_cvp.png"), p, width = 12,
           height = max(6, 3 * length(unique(subjects$analyte))), dpi = 200)
  }
  if (nrow(sensitivity)) {
    sensitivity$scenario <- factor(sensitivity$scenario, levels = rev(unique(sensitivity$scenario)))
    p <- ggplot(sensitivity, aes(x = cvi_median_pct, y = scenario)) +
      geom_segment(aes(x = cvi_q025_pct, xend = cvi_q975_pct, yend = scenario), colour = "#547795") +
      geom_point(colour = "#153854") + facet_wrap(~analyte + group, scales = "free_x", ncol = 2) +
      labs(x = "Population CVI (%)", y = NULL, title = "Prior sensitivity",
           subtitle = "Identical observations in every scenario; posterior median and 95% CrI")
    ggsave(file.path(out, "prior_sensitivity.png"), p, width = 12,
           height = max(6, 3 * length(unique(sensitivity$analyte))), dpi = 200)
  }
}
