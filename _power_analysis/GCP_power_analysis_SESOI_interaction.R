## Classical (SIMR) power for Contrast × Microsaccade interaction with residual-SD sensitivity.
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

  sesoi_beta_interaction <- -0.10
  true_beta_interaction <- 1.5 * sesoi_beta_interaction
  pilot_beta_interaction <- 0.0100920695330325

  true_beta_contrast <- 0.00225168622867046
  true_beta_microsaccade <- -0.00719377561756273
  outcome_mean <- 0.370032495163686
  microsaccade_mean <- -9.81325075598174

  baseline_random_intercept_sd <- 0.332724620668064
  # Pilot N=10 weakly identifies slope variances; apply conservative non-zero floors.
  baseline_random_slope_contrast_sd <- 0.04
  baseline_random_slope_microsaccade_sd <- 0.04

  baseline_ms_subject_sd <- 27.4543840173968
  baseline_ms_residual_sd <- 86.4566573636538
  true_ms_effect_contrast <- -5.58866992826615

  residual_sd_levels <- c(0.1932, 0.2976, 0.4230)
  baseline_residual_sd <- residual_sd_levels[2L]
  residual_multipliers <- residual_sd_levels / baseline_residual_sd
  log_msg <- function(...) {
    cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = "")
    flush.console()
  }

  output_prefix <- "GCP_power_analysis_SESOI_interaction"
  figure_res_dpi <- 600L
  plot_title <- sprintf(
    "SESOI interaction classical power (SESOI=%s; true(DGP)=%s; pilot=%s)",
    format(sesoi_beta_interaction, digits = 3),
    format(true_beta_interaction, digits = 3),
    format(pilot_beta_interaction, digits = 3)
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

  simulate_microsaccade <- function(dat) {
    n_sub <- nlevels(dat$Subject)
    ms_intercept_subj <- rnorm(n_sub, mean = 0, sd = baseline_ms_subject_sd)
    ms_raw <- microsaccade_mean +
      true_ms_effect_contrast * dat$contrast_num_c +
      ms_intercept_subj[dat$Subject] +
      rnorm(nrow(dat), mean = 0, sd = baseline_ms_residual_sd)
    as.numeric(scale(ms_raw, center = TRUE, scale = TRUE))
  }

  simulate_response <- function(dat, residual_sd) {
    n_sub <- nlevels(dat$Subject)
    u0 <- rnorm(n_sub, mean = 0, sd = baseline_random_intercept_sd)
    u1 <- rnorm(n_sub, mean = 0, sd = baseline_random_slope_contrast_sd)
    u2 <- rnorm(n_sub, mean = 0, sd = baseline_random_slope_microsaccade_sd)
    x <- dat$contrast_num_c
    m <- dat$microsaccade_c
    y_mu <- outcome_mean +
      u0[dat$Subject] +
      (true_beta_contrast + u1[dat$Subject]) * x +
      (true_beta_microsaccade + u2[dat$Subject]) * m +
      true_beta_interaction * x * m
    y_mu + rnorm(nrow(dat), mean = 0, sd = residual_sd)
  }

  set.seed(seed)
  log_msg("Building simulation template and base model.")
  template_dat <- make_template(max(subject_breaks))
  template_dat$microsaccade_c <- simulate_microsaccade(template_dat)
  template_dat$gamma_power <- simulate_response(template_dat, residual_sd = baseline_residual_sd)

  base_model <- lmer(
    gamma_power ~ contrast_num_c * microsaccade_c + (1 + contrast_num_c + microsaccade_c | Subject),
    data = template_dat,
    REML = FALSE,
    control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
  )

  fixef(base_model)["(Intercept)"] <- outcome_mean
  fixef(base_model)["contrast_num_c"] <- true_beta_contrast
  fixef(base_model)["microsaccade_c"] <- true_beta_microsaccade
  fixef(base_model)["contrast_num_c:microsaccade_c"] <- true_beta_interaction
  VarCorr(base_model) <- list(
    Subject = matrix(
      c(
        baseline_random_intercept_sd^2, 0, 0,
        0, baseline_random_slope_contrast_sd^2, 0,
        0, 0, baseline_random_slope_microsaccade_sd^2
      ),
      nrow = 3,
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
    log_msg(sprintf("Running scenario %d/%d | residual SD=%s | nsim=%d", i, nrow(scenario_df), sc$residual_sd_label, nsim))
    model_sc <- base_model
    sigma(model_sc) <- sc$residual_sd
    pc <- suppressWarnings(powerCurve(
      model_sc,
      test = fixed("contrast_num_c:microsaccade_c"),
      along = "Subject",
      breaks = subject_breaks,
      nsim = nsim,
      progress = TRUE
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
  log_msg("Saved curve CSV.")

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
  log_msg("Saved summary CSV.")

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
  log_msg("Saved line plot PNG.")

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
  log_msg("Saved heatmap PNG.")

  baseline_model <- base_model
  sigma(baseline_model) <- baseline_residual_sd
  baseline_pc <- suppressWarnings(powerCurve(
    baseline_model,
    test = fixed("contrast_num_c:microsaccade_c"),
    along = "Subject",
    breaks = subject_breaks,
    nsim = nsim,
    progress = TRUE
  ))
  png(
    file = file.path(figure_output_dir, paste0(output_prefix, "_simr_powercurve.png")),
    width = 10, height = 6.5, units = "in", res = figure_res_dpi
  )
  plot(baseline_pc)
  dev.off()
  log_msg("Saved SIMR default powerCurve PNG.")

  power_df
}

result <- runSESOI()
