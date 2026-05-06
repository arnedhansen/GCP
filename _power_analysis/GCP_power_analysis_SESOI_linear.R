## Classical (SIMR) power for linear contrast effect with residual-SD sensitivity.
## Outputs: CSV, line plot, heatmap, and additional default SIMR powerCurve plot.

ensure_packages <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    install.packages(missing, repos = "https://cloud.r-project.org")
  }
}

ensure_packages(c("lme4", "simr", "ggplot2", "scales"))
suppressPackageStartupMessages({
  library(lme4)
  library(simr)
  library(ggplot2)
  library(scales)
})

resolve_gcp_root <- function() {
  if (.Platform$OS.type == "windows") {
    return("W:/Students/Arne/GCP")
  }
  "/Volumes/g_psyplafor_methlab$/Students/Arne/GCP"
}

extract_powercurve_df <- function(pc, n_subjects_vec) {
  do.call(rbind, lapply(seq_along(pc$ps), function(i) {
    z <- pc$ps[[i]]
    ci <- suppressWarnings(confint(z))
    data.frame(
      n_subjects = n_subjects_vec[i],
      power = as.numeric(z$x / z$n),
      lower = as.numeric(ci[1, 1]),
      upper = as.numeric(ci[1, 2]),
      nsim = as.integer(z$n),
      stringsAsFactors = FALSE
    )
  }))
}

runSESOI <- function() {
  seed <- 123
  nsim <- 100
  strict_power_target <- 0.90
  subject_breaks <- seq(10, 100, by = 10)
  contrast_levels <- c("25", "50", "75", "100")
  trials_per_condition <- 160

  sesoi_beta <- 0.10
  true_beta <- 1.5 * sesoi_beta
  pilot_true_beta_contrast_num_c <- 0.247639934007685

  outcome_mean <- 53.7513302255983
  baseline_random_intercept_sd <- 2.95596738123282
  pilot_random_slope_sd <- 0.00138655105795991
  # Pilot N=10 yields near-zero slope variance; apply conservative non-zero floor.
  baseline_random_slope_sd <- max(pilot_random_slope_sd, 0.05)

  residual_sd_levels <- c(0.4645, 0.6935, 0.8923)
  baseline_residual_sd <- residual_sd_levels[2L]
  residual_multipliers <- residual_sd_levels / baseline_residual_sd

  output_prefix <- "GCP_power_analysis_SESOI_linear"
  figure_res_dpi <- 600L
  plot_title <- sprintf(
    "SESOI linear classical power (SESOI=%s; true(DGP)=%s; pilot c=%s)",
    format(sesoi_beta, digits = 3),
    format(true_beta, digits = 3),
    format(pilot_true_beta_contrast_num_c, digits = 3)
  )

  gcp_root <- resolve_gcp_root()
  figure_output_dir <- file.path(gcp_root, "figures", "power_analysis")
  data_output_dir <- file.path(gcp_root, "data", "power_analysis")
  dir.create(figure_output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(data_output_dir, recursive = TRUE, showWarnings = FALSE)

  make_template <- function(n_subjects) {
    Subject <- factor(rep(seq_len(n_subjects), each = length(contrast_levels) * trials_per_condition))
    contrast <- factor(
      rep(rep(contrast_levels, each = trials_per_condition), times = n_subjects),
      levels = contrast_levels,
      ordered = TRUE
    )
    contrast_num <- as.numeric(as.character(contrast))
    contrast_num_c <- as.numeric(scale(contrast_num, center = TRUE, scale = TRUE))
    data.frame(Subject = Subject, contrast_num_c = contrast_num_c)
  }

  simulate_response <- function(dat, beta, ri_sd, rs_sd, residual_sd, mu) {
    n_sub <- nlevels(dat$Subject)
    u0 <- rnorm(n_sub, mean = 0, sd = ri_sd)
    u1 <- rnorm(n_sub, mean = 0, sd = rs_sd)
    x <- dat$contrast_num_c
    y_mu <- mu + u0[dat$Subject] + (beta + u1[dat$Subject]) * x
    y_mu + rnorm(nrow(dat), mean = 0, sd = residual_sd)
  }

  set.seed(seed)
  template_dat <- make_template(max(subject_breaks))
  template_dat$gamma_power <- simulate_response(
    template_dat,
    beta = true_beta,
    ri_sd = baseline_random_intercept_sd,
    rs_sd = baseline_random_slope_sd,
    residual_sd = baseline_residual_sd,
    mu = outcome_mean
  )

  base_model <- lmer(
    gamma_power ~ contrast_num_c + (1 + contrast_num_c | Subject),
    data = template_dat,
    REML = FALSE,
    control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
  )

  fixef(base_model)["(Intercept)"] <- outcome_mean
  fixef(base_model)["contrast_num_c"] <- true_beta
  VarCorr(base_model) <- list(
    Subject = matrix(
      c(
        baseline_random_intercept_sd^2, 0,
        0, baseline_random_slope_sd^2
      ),
      nrow = 2,
      byrow = TRUE
    )
  )

  scenario_df <- data.frame(
    scenario_label = sprintf("res_%.2f", residual_multipliers),
    residual_sd = residual_sd_levels,
    residual_sd_label = sprintf("%.2f", residual_sd_levels),
    stringsAsFactors = FALSE
  )

  power_rows <- lapply(seq_len(nrow(scenario_df)), function(i) {
    sc <- scenario_df[i, , drop = FALSE]
    model_sc <- base_model
    sigma(model_sc) <- sc$residual_sd
    pc <- suppressWarnings(powerCurve(
      model_sc,
      test = fixed("contrast_num_c"),
      along = "Subject",
      breaks = subject_breaks,
      nsim = nsim,
      progress = FALSE
    ))
    out <- extract_powercurve_df(pc, subject_breaks)
    out$scenario_label <- sc$scenario_label
    out$residual_sd <- sc$residual_sd
    out$residual_sd_label <- sc$residual_sd_label
    out
  })
  power_df <- do.call(rbind, power_rows)

  power_df$power_unconditional <- power_df$power
  power_df$lower_unconditional <- power_df$lower
  power_df$upper_unconditional <- power_df$upper

  write.csv(power_df, file.path(data_output_dir, paste0(output_prefix, "_curve.csv")), row.names = FALSE)

  summary_rows <- do.call(rbind, lapply(split(power_df, power_df$scenario_label), function(df) {
    df <- df[order(df$n_subjects), , drop = FALSE]
    hit <- df[df$power_unconditional >= strict_power_target, "n_subjects", drop = TRUE]
    data.frame(
      scenario_label = as.character(df$scenario_label[1]),
      residual_sd = as.numeric(df$residual_sd[1]),
      N_min_for_90 = if (length(hit) > 0) min(hit) else NA_integer_,
      power_at_max_N = df$power_unconditional[which.max(df$n_subjects)],
      stringsAsFactors = FALSE
    )
  }))
  write.csv(summary_rows, file.path(data_output_dir, paste0(output_prefix, "_summary.csv")), row.names = FALSE)

  curve_plot <- ggplot(
    power_df,
    aes(x = .data$n_subjects, y = .data$power_unconditional, color = .data$residual_sd_label, group = .data$residual_sd_label)
  ) +
    geom_errorbar(aes(ymin = .data$lower_unconditional, ymax = .data$upper_unconditional), linewidth = 0.6, width = 1.6, alpha = 0.70) +
    geom_line(linewidth = 0.9, linetype = "dotted") +
    geom_point(size = 2.8) +
    geom_hline(yintercept = strict_power_target, linetype = "dashed", color = "grey45", linewidth = 0.7) +
    scale_color_manual(
      values = stats::setNames(c("#1B9E77", "#D95F02", "#7570B3"), sprintf("%.2f", residual_sd_levels)),
      breaks = sprintf("%.2f", residual_sd_levels),
      labels = sprintf("%.2f", residual_sd_levels)
    ) +
    scale_x_continuous(breaks = sort(unique(subject_breaks))) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1), breaks = c(0, 0.25, 0.50, 0.75, 0.90, 1)) +
    labs(x = "Subjects", y = "Classical power", color = "Residuals SD", title = plot_title) +
    theme_classic(base_size = 18, base_family = "Arial") +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      panel.grid = element_blank(),
      axis.line = element_line(color = "black", linewidth = 1.2),
      axis.ticks = element_line(color = "black", linewidth = 1.0),
      axis.ticks.length = grid::unit(0.25, "cm"),
      axis.title = element_text(size = 20),
      axis.text = element_text(size = 16, color = "black"),
      legend.position = "right",
      legend.title = element_text(size = 15),
      legend.text = element_text(size = 13),
      plot.title = element_text(size = 24, face = "plain", hjust = 0.5)
    )

  png(
    file = file.path(figure_output_dir, paste0(output_prefix, ".png")),
    width = 10, height = 6.5, units = "in", res = figure_res_dpi
  )
  print(curve_plot)
  dev.off()

  heatmap_plot <- ggplot(power_df, aes(x = factor(.data$n_subjects), y = .data$residual_sd_label, fill = .data$power_unconditional)) +
    geom_tile(color = "white", linewidth = 1.1) +
    geom_text(aes(label = sprintf("%.2f", .data$power_unconditional)), color = "white", size = 3.5, fontface = "bold", family = "Arial") +
    scale_fill_gradientn(
      colours = c("#8E0F1F", "#F46D43", "#FEE08B", "#66BD63", "#0B6E3A"),
      values = c(0.00, 0.60, 0.79, 0.80, 1.00),
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.2),
      labels = sprintf("%.2f", seq(0, 1, by = 0.2))
    ) +
    coord_fixed() +
    labs(x = "Subjects", y = "Residuals SD", fill = "Classical power", title = plot_title) +
    theme_minimal(base_size = 16, base_family = "Arial") +
    theme(
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      axis.title = element_text(size = 13),
      axis.text = element_text(size = 11),
      legend.position = "right",
      legend.title = element_text(size = 11),
      legend.text = element_text(size = 10),
      plot.title = element_text(size = 20, face = "plain", hjust = 0.5)
    )

  png(
    file = file.path(figure_output_dir, paste0(output_prefix, "_heatmap.png")),
    width = 11, height = 6.5, units = "in", res = figure_res_dpi
  )
  print(heatmap_plot)
  dev.off()

  baseline_model <- base_model
  sigma(baseline_model) <- baseline_residual_sd
  baseline_pc <- suppressWarnings(powerCurve(
    baseline_model,
    test = fixed("contrast_num_c"),
    along = "Subject",
    breaks = subject_breaks,
    nsim = nsim,
    progress = FALSE
  ))
  png(
    file = file.path(figure_output_dir, paste0(output_prefix, "_simr_powercurve.png")),
    width = 10, height = 6.5, units = "in", res = figure_res_dpi
  )
  plot(baseline_pc)
  dev.off()

  power_df
}

result <- runSESOI()
