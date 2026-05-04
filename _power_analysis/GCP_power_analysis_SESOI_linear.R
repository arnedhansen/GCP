## Simulates SESOI-based power with mixed models and saves a CSV power curve plus PNG figure

# install.packages("lme4", repos = "https://cloud.r-project.org")
# install.packages("ggplot2", repos = "https://cloud.r-project.org")
# install.packages("scales", repos = "https://cloud.r-project.org")

# PeakFrequency pilot (GCP_power_analysis_pilot_stats.R): fixed effects from pilot_mixed_model_fixed_effects_*;
# VarCorr from pilot_mixed_model_random_effects_* after refit with (1|Subject) because (1+contrast|Subject) was singular.
# Random slope SD is therefore not identified in that export (set to 0 in the DGP below).
# Residual SD scenarios: bootstrap q2.5/q50/q97.5 of sigma (subject resampling, same lmer formula as final pilot model);
#   see pilot_subject_level_residual_sigma_bootstrap.csv (sync after rerunning pilot_stats).
# Source: .../data/pilot_stats/

suppressPackageStartupMessages(library(lme4))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(scales))

resolve_gcp_root <- function() {
  if (.Platform$OS.type == "windows") {
    return("W:/Students/Arne/GCP")
  }
  "/Volumes/g_psyplafor_methlab$/Students/Arne/GCP"
}

runSESOI <- function() {
  # Simulation settings
  seed <- 123
  sesoi_decision_alpha <- 0.05
  nsim <- 100
  strict_power_target <- 0.90
  subject_breaks <- seq(10, 100, by = 10)
  contrast_levels <- c("25", "50", "75", "100")
  trials_per_condition <- 160
  sesoi_beta <- 0.10
  true_beta <- 0.247639934007685 ## pilot contrast_num_c (PeakFrequency), simplified lmer
  outcome_mean <- 53.7513302255983 ## power_parameter_manifest.csv (PeakFrequency)
  baseline_random_intercept_sd <- 2.95596738123282 ## pilot VarCorr Subject (Intercept); simplified lmer
  baseline_random_slope_sd <- 0 ## no slope variance in pilot RI-only refit
  residual_sd_levels <- c(0.507149145838165, 0.719094938085036, 0.893798456396689) ## pilot_subject_level_residual_sigma_bootstrap.csv (PeakFrequency)
  baseline_residual_sd <- residual_sd_levels[2L] ## bootstrap median sigma (scenarios are q2.5/median/q97.5)
  ri_multiplier_fixed <- 1.00
  residual_multipliers <- residual_sd_levels / baseline_residual_sd
  trial_missingness_rate <- 0.20
  subject_dropout_rate <- 0.10
  parallel_workers <- 8
  parallel_round_chunk_nsim <- 25
  sim_col <- "gamma_power"
  model_formula_full <- gamma_power ~ contrast_num_c + (1 + contrast_num_c | Subject)
  output_prefix <- "GCP_power_analysis_SESOI_linear"
  figure_res_dpi <- 600L
  plot_title <- "Power Analysis: SESOI Linear Slope"
  power_label <- "Power"
  # Darker red-green endpoints with green onset at >= 0.8 in heatmaps
  heat_low_color <- "#8E0F1F"
  heat_high_color <- "#0B6E3A"

  gcp_root <- resolve_gcp_root()
  figure_output_dir <- file.path(gcp_root, "figures", "power_analysis")
  data_output_dir <- file.path(gcp_root, "data", "power_analysis")
  curve_csv_path <- file.path(data_output_dir, paste0(output_prefix, "_curve.csv"))
  set.seed(seed)
  dir.create(figure_output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(data_output_dir, recursive = TRUE, showWarnings = FALSE)

  # Build compact scenario grid varying only residual SD
  make_scenarios <- function() {
    scenario_df <- data.frame(
      residual_multiplier = residual_multipliers,
      stringsAsFactors = FALSE
    )
    scenario_df$scenario_label <- sprintf(
      "res_%.2f",
      scenario_df$residual_multiplier
    )
    scenario_df$varied_component <- "residual_only"
    scenario_df$multiplier <- scenario_df$residual_multiplier
    scenario_df$random_intercept_sd_value <- baseline_random_intercept_sd * ri_multiplier_fixed
    scenario_df$random_slope_sd_value <- baseline_random_slope_sd
    scenario_df$residual_sd_value <- baseline_residual_sd * scenario_df$residual_multiplier
    scenario_df
  }

  estimate_power_chunk <- function(
      chunk_nsim,
      n_subjects,
      trials_per_condition,
      sim_col,
      model_formula_full,
      sesoi_decision_alpha,
      beta_raw,
      true_beta,
      outcome_mean,
      random_intercept_sd,
      random_slope_sd,
      residual_sd,
      trial_missingness_rate,
      subject_dropout_rate) {
    # Count SESOI-decision outcomes and diagnostic metrics
    fit_attempts <- 0
    valid_fits <- 0
    sesoi_successes <- 0
    sign_errors <- 0
    type_m_sum <- 0
    type_m_count <- 0
    singular_count <- 0
    convergence_warn_count <- 0
    for (i in seq_len(chunk_nsim)) {
      # Create balanced trial table before missingness and dropout
      Subject <- factor(rep(seq_len(n_subjects), each = length(contrast_levels) * trials_per_condition))
      contrast <- factor(
        rep(rep(contrast_levels, each = trials_per_condition), times = n_subjects),
        levels = contrast_levels,
        ordered = TRUE
      )
      contrast_num <- as.numeric(as.character(contrast))
      # Standardize contrast predictor
      contrast_num_c <- as.numeric(scale(contrast_num, center = TRUE, scale = TRUE))
      dat <- data.frame(Subject = Subject, contrast_num_c = contrast_num_c)

      # Generate outcome with fixed SESOI, subject random effects, and residual noise
      n_subject_levels <- nlevels(dat$Subject)
      random_intercepts <- rnorm(n_subject_levels, mean = 0, sd = random_intercept_sd)
      random_slopes <- rnorm(n_subject_levels, mean = 0, sd = random_slope_sd)
      x <- dat$contrast_num_c
      mu <- outcome_mean + random_intercepts[dat$Subject] + random_slopes[dat$Subject] * x + true_beta * x
      dat[[sim_col]] <- mu + rnorm(nrow(dat), mean = 0, sd = residual_sd)

      # Apply trial-level missingness
      keep_trial <- stats::runif(nrow(dat)) > trial_missingness_rate
      dat <- dat[keep_trial, , drop = FALSE]
      if (nrow(dat) == 0) {
        next
      }

      # Apply subject-level dropout after trial loss
      subject_levels <- levels(droplevels(dat$Subject))
      n_drop_subjects <- floor(length(subject_levels) * subject_dropout_rate)
      if (n_drop_subjects > 0) {
        dropped_subjects <- sample(subject_levels, size = n_drop_subjects, replace = FALSE)
        dat <- dat[!(dat$Subject %in% dropped_subjects), , drop = FALSE]
      }
      dat <- droplevels(dat)
      if (nlevels(dat$Subject) < 3 || nrow(dat) == 0) {
        next
      }

      # Fit nested models and compare by likelihood-ratio test
      fit_full <- tryCatch(
        suppressMessages(lmer(model_formula_full, data = dat, REML = FALSE)),
        error = function(e) NULL
      )
      if (is.null(fit_full)) {
        next
      }
      coef_tbl <- tryCatch(summary(fit_full)$coefficients, error = function(e) NULL)
      if (is.null(coef_tbl) || !"contrast_num_c" %in% rownames(coef_tbl)) {
        next
      }
      fit_attempts <- fit_attempts + 1
      is_singular_fit <- isTRUE(lme4::isSingular(fit_full, tol = 1e-4))
      singular_count <- singular_count + as.integer(is_singular_fit)
      if (is_singular_fit) {
        next
      }
      valid_fits <- valid_fits + 1
      conv_msgs <- fit_full@optinfo$conv$lme4$messages
      convergence_warn_count <- convergence_warn_count + as.integer(length(conv_msgs) > 0)

      est <- as.numeric(coef_tbl["contrast_num_c", "Estimate"])
      est_se <- as.numeric(coef_tbl["contrast_num_c", "Std. Error"])
      z_one_sided <- stats::qnorm(1 - sesoi_decision_alpha)
      ci_low_one_sided <- est - z_one_sided * est_se
      ci_high_one_sided <- est + z_one_sided * est_se
      sesoi_success <- if (beta_raw >= 0) ci_low_one_sided > beta_raw else ci_high_one_sided < beta_raw
      sesoi_successes <- sesoi_successes + as.integer(sesoi_success)
      sign_errors <- sign_errors + as.integer(sign(est) != sign(true_beta))
      if (isTRUE(sesoi_success)) {
        type_m_sum <- type_m_sum + abs(est / true_beta)
        type_m_count <- type_m_count + 1
      }
    }
    list(
      chunk_nsim = as.integer(chunk_nsim),
      fit_attempts = as.integer(fit_attempts),
      valid_fits = as.integer(valid_fits),
      sesoi_successes = as.integer(sesoi_successes),
      sign_errors = as.integer(sign_errors),
      type_m_sum = as.numeric(type_m_sum),
      type_m_count = as.integer(type_m_count),
      singular_count = as.integer(singular_count),
      convergence_warn_count = as.integer(convergence_warn_count)
    )
  }

  cl <- parallel::makeCluster(as.integer(parallel_workers), type = "PSOCK")
  on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::clusterEvalQ(cl, suppressPackageStartupMessages(library(lme4)))
  parallel::clusterSetRNGStream(cl, iseed = seed)
  parallel::clusterExport(cl, varlist = c("estimate_power_chunk", "contrast_levels"), envir = environment())

  scenario_df <- make_scenarios()
  scenario_order <- scenario_df$scenario_label

  scenario_rows <- lapply(seq_len(nrow(scenario_df)), function(idx) {
    scenario <- scenario_df[idx, , drop = FALSE]
    n_rows <- lapply(subject_breaks, function(n_subjects) {
      cat(sprintf("[%s] START scenario=%s N=%d", format(Sys.time(), "%H:%M:%S"), scenario$scenario_label, n_subjects))
      flush.console()

      total_nsim <- 0
      total_fit_attempts <- 0
      total_valid_fits <- 0
      total_sesoi_successes <- 0
      total_sign_errors <- 0
      total_type_m_sum <- 0
      total_type_m_count <- 0
      total_singular <- 0
      total_convergence_warn <- 0
      # Run chunks in parallel until nsim is reached
      while (total_nsim < nsim) {
        remaining <- nsim - total_nsim
        workers_this_round <- min(length(cl), remaining)
        chunk_sizes <- rep(parallel_round_chunk_nsim, workers_this_round)
        overflow <- (workers_this_round * parallel_round_chunk_nsim) - remaining
        if (overflow > 0) {
          for (k in seq_len(overflow)) {
            idx2 <- ((k - 1) %% workers_this_round) + 1
            chunk_sizes[idx2] <- chunk_sizes[idx2] - 1
          }
        }
        chunk_sizes <- chunk_sizes[chunk_sizes > 0]

        chunk_results <- parallel::parLapply(
          cl = cl,
          X = chunk_sizes,
          fun = estimate_power_chunk,
          n_subjects = n_subjects,
          trials_per_condition = trials_per_condition,
          sim_col = sim_col,
          model_formula_full = model_formula_full,
          sesoi_decision_alpha = sesoi_decision_alpha,
          beta_raw = sesoi_beta,
          true_beta = true_beta,
          outcome_mean = outcome_mean,
          random_intercept_sd = scenario$random_intercept_sd_value,
          random_slope_sd = scenario$random_slope_sd_value,
          residual_sd = scenario$residual_sd_value,
          trial_missingness_rate = trial_missingness_rate,
          subject_dropout_rate = subject_dropout_rate
        )

        round_n <- sum(vapply(chunk_results, function(x) as.integer(x$chunk_nsim), integer(1)))
        round_fit_attempts <- sum(vapply(chunk_results, function(x) as.integer(x$fit_attempts), integer(1)))
        round_valid_fits <- sum(vapply(chunk_results, function(x) as.integer(x$valid_fits), integer(1)))
        round_sesoi_successes <- sum(vapply(chunk_results, function(x) as.integer(x$sesoi_successes), integer(1)))
        round_sign_errors <- sum(vapply(chunk_results, function(x) as.integer(x$sign_errors), integer(1)))
        round_type_m_sum <- sum(vapply(chunk_results, function(x) as.numeric(x$type_m_sum), numeric(1)))
        round_type_m_count <- sum(vapply(chunk_results, function(x) as.integer(x$type_m_count), integer(1)))
        round_singular <- sum(vapply(chunk_results, function(x) as.integer(x$singular_count), integer(1)))
        round_convergence_warn <- sum(vapply(chunk_results, function(x) as.integer(x$convergence_warn_count), integer(1)))
        total_nsim <- total_nsim + round_n
        total_fit_attempts <- total_fit_attempts + round_fit_attempts
        total_valid_fits <- total_valid_fits + round_valid_fits
        total_sesoi_successes <- total_sesoi_successes + round_sesoi_successes
        total_sign_errors <- total_sign_errors + round_sign_errors
        total_type_m_sum <- total_type_m_sum + round_type_m_sum
        total_type_m_count <- total_type_m_count + round_type_m_count
        total_singular <- total_singular + round_singular
        total_convergence_warn <- total_convergence_warn + round_convergence_warn
      }

      power_conditional <- if (total_valid_fits > 0) total_sesoi_successes / total_valid_fits else NA_real_
      power_unconditional <- if (total_nsim > 0) total_sesoi_successes / total_nsim else NA_real_
      valid_fit_rate <- if (total_nsim > 0) total_valid_fits / total_nsim else NA_real_
      fit_failure_rate <- if (total_nsim > 0) 1 - valid_fit_rate else NA_real_
      se_cond <- if (total_valid_fits > 0) {
        sqrt(power_conditional * (1 - power_conditional) / total_valid_fits)
      } else {
        NA_real_
      }
      se_uncond <- if (total_nsim > 0) {
        sqrt(power_unconditional * (1 - power_unconditional) / total_nsim)
      } else {
        NA_real_
      }
      cat(sprintf(
        " | cond=%s uncond=%s valid_fit=%.3f fit_fail=%.3f\n",
        if (is.finite(power_conditional)) sub("^0", "", sprintf("%.3f", power_conditional)) else "NA",
        if (is.finite(power_unconditional)) sub("^0", "", sprintf("%.3f", power_unconditional)) else "NA",
        valid_fit_rate,
        fit_failure_rate
      ))
      flush.console()
      out <- data.frame(
        scenario_label = scenario$scenario_label,
        varied_component = scenario$varied_component,
        multiplier = scenario$multiplier,
        random_intercept_sd = scenario$random_intercept_sd_value,
        random_slope_sd = scenario$random_slope_sd_value,
        residual_sd = scenario$residual_sd_value,
        n_subjects = n_subjects,
        power = power_conditional,
        power_conditional = power_conditional,
        power_unconditional = power_unconditional,
        lower = if (is.finite(se_cond)) pmax(0, power_conditional - 1.96 * se_cond) else NA_real_,
        upper = if (is.finite(se_cond)) pmin(1, power_conditional + 1.96 * se_cond) else NA_real_,
        lower_unconditional = if (is.finite(se_uncond)) {
          pmax(0, power_unconditional - 1.96 * se_uncond)
        } else {
          NA_real_
        },
        upper_unconditional = if (is.finite(se_uncond)) {
          pmin(1, power_unconditional + 1.96 * se_uncond)
        } else {
          NA_real_
        },
        nsim = total_nsim,
        valid_fits = total_valid_fits,
        valid_fit_rate = valid_fit_rate,
        fit_failure_rate = fit_failure_rate,
        fit_success_rate = valid_fit_rate,
        type_s = if (total_valid_fits > 0) total_sign_errors / total_valid_fits else NA_real_,
        type_m = if (total_type_m_count > 0) total_type_m_sum / total_type_m_count else NA_real_,
        singular_rate = if (total_fit_attempts > 0) total_singular / total_fit_attempts else NA_real_,
        fit_attempt_rate = if (total_nsim > 0) total_fit_attempts / total_nsim else NA_real_,
        convergence_warn_rate = if (total_valid_fits > 0) total_convergence_warn / total_valid_fits else NA_real_,
        stringsAsFactors = FALSE
      )
      out
    })
    do.call(rbind, n_rows)
  })

  power_df <- do.call(rbind, scenario_rows)
  power_df$scenario_label <- factor(power_df$scenario_label, levels = scenario_order, ordered = TRUE)
  power_df$meets_target_90 <- power_df$power_conditional >= strict_power_target
  write.csv(power_df, curve_csv_path, row.names = FALSE)

  power_df$residual_sd_label <- factor(
    sprintf("%.2f", power_df$residual_sd),
    levels = sprintf("%.2f", residual_sd_levels)
  )
  approx_equal <- function(x, y, tol = 1e-8) {
    abs(x - y) < tol
  }
  stopifnot(length(unique(power_df$scenario_label)) == length(residual_multipliers))
  stopifnot(all(approx_equal(power_df$random_intercept_sd, baseline_random_intercept_sd * ri_multiplier_fixed)))
  stopifnot(all(approx_equal(power_df$random_slope_sd, baseline_random_slope_sd)))

  # Summarize minimum sample size for target power
  summary_rows <- do.call(rbind, lapply(split(power_df, power_df$scenario_label), function(df) {
    df <- df[order(df$n_subjects), , drop = FALSE]
    hit <- df[df$meets_target_90, "n_subjects", drop = TRUE]
    max_idx <- which.max(df$n_subjects)
    data.frame(
      scenario_label = as.character(df$scenario_label[1]),
      varied_component = as.character(df$varied_component[1]),
      multiplier = as.numeric(df$multiplier[1]),
      N_min_for_90_conditional = if (length(hit) > 0) min(hit) else NA_integer_,
      power_conditional_at_max_N = df$power_conditional[max_idx],
      power_unconditional_at_max_N = df$power_unconditional[max_idx],
      valid_fit_rate_at_max_N = df$valid_fit_rate[max_idx],
      fit_failure_rate_at_max_N = df$fit_failure_rate[max_idx],
      type_s_at_max_N = df$type_s[max_idx],
      type_m_at_max_N = df$type_m[max_idx],
      singular_rate_at_max_N = df$singular_rate[max_idx],
      convergence_warn_rate_at_max_N = df$convergence_warn_rate[max_idx],
      monotonic_non_decreasing_conditional = all(diff(df$power_conditional) >= -0.02),
      stringsAsFactors = FALSE
    )
  }))
  write.csv(summary_rows, file.path(data_output_dir, paste0(output_prefix, "_summary.csv")), row.names = FALSE)
  curve_plot <- ggplot(power_df, aes(x = .data$n_subjects, y = .data$power_conditional, color = .data$residual_sd_label, group = .data$residual_sd_label)) +
    geom_errorbar(aes(ymin = .data$lower, ymax = .data$upper), linewidth = 0.6, width = 1.6, alpha = 0.70) +
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
    labs(x = "Subjects", y = power_label, color = "Residuals SD", title = plot_title) +
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
  png(file = file.path(figure_output_dir, paste0(output_prefix, ".png")), width = 2200, height = 1400, res = figure_res_dpi)
  print(curve_plot)
  dev.off()

  # Plot single-block heatmap with residual rows and subject columns
  heatmap_plot <- ggplot(power_df, aes(x = factor(.data$n_subjects), y = .data$residual_sd_label, fill = .data$power_conditional)) +
    geom_tile(color = "white", linewidth = 1.1) +
    geom_text(aes(label = sprintf("%.2f", .data$power_conditional)), color = "white", size = 3.5, fontface = "bold", family = "Arial") +
    scale_fill_gradientn(
      colours = c(heat_low_color, "#F46D43", "#FEE08B", "#66BD63", heat_high_color),
      values = c(0.00, 0.60, 0.79, 0.80, 1.00),
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.2),
      labels = sprintf("%.2f", seq(0, 1, by = 0.2))
    ) +
    coord_fixed() +
    labs(
      x = "Subjects",
      y = "Residuals SD",
      fill = power_label,
      title = plot_title
    ) +
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
      plot.title = element_text(size = 20, face = "plain", hjust = 0.5),
      panel.spacing = grid::unit(0.9, "lines")
    )
  png(file = file.path(figure_output_dir, paste0(output_prefix, "_heatmap.png")), width = 2500, height = 1400, res = figure_res_dpi)
  print(heatmap_plot)
  dev.off()

  power_df
}

result <- runSESOI()