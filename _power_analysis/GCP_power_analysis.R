## GCP simulation-based power analysis (Stage 1).
## Models: H5 linear gamma frequency, H6 quadratic gamma power, H7 contrast x microsaccade interaction.
## Detection power only (two-sided alpha = 0.05). Recommended N from pessimistic scenarios.

ensure_packages <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    install.packages(missing, repos = "https://cloud.r-project.org")
  }
}

ensure_packages(c("lme4", "lmerTest", "ggplot2", "scales", "simr"))
suppressPackageStartupMessages({
  library(lme4)
  library(lmerTest)
  library(ggplot2)
  library(scales)
  library(simr)
})

script_dir <- {
  cmd <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd, value = TRUE)
  if (length(file_arg) > 0) {
    dirname(sub("^--file=", "", file_arg[1]))
  } else {
    file.path(getwd(), "_power_analysis")
  }
}

source(file.path(script_dir, "GCP_power_analysis_helpers.R"))
source(file.path(script_dir, "GCP_power_analysis_engine.R"))

run_power_analysis <- function(
    nsim_subject = 1000L,
    nsim_interaction = 500L,
    subject_breaks = seq(20, 100, by = 10),
    parallel_workers = 8L,
    power_target = 0.90,
    plot_only = FALSE) {
  gcp_root <- resolve_gcp_root()
  data_dir <- file.path(gcp_root, "data", "power_analysis")
  dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

  partition_freq <- load_variance_partition_summary(gcp_root, "PeakFrequency", "linear")
  partition_power <- load_variance_partition_summary(gcp_root, "PeakAmplitude", "quadratic")

  beta_freq <- load_conservative_simulation_beta(
    gcp_root, "PeakFrequency", "contrast_num_c", direction = "positive", pilot_fallback = 0.05
  )
  beta_quad <- load_conservative_simulation_beta(
    gcp_root, "PeakAmplitude", "contrast_num_c2", direction = "negative", pilot_fallback = -0.04
  )
  beta_linear_nuisance <- load_pilot_fixed_effect(
    gcp_root, "PeakAmplitude", "contrast_num_c", fallback = 0.077
  )

  message(sprintf("Simulation betas (conservative, not Table 1): H5 beta=%.4f | H6 beta=%.4f", beta_freq, beta_quad))

  message("=== H5: Linear gamma peak frequency (contrast_num_c) ===")
  message("H5 | residual sensitivity (3 residual-SD scenarios)")
  res_freq <- run_subject_level_power(
    config = list(
      label = "Gamma peak frequency (linear contrast)",
      hypothesis = "H5",
      pilot_outcome = "PeakFrequency",
      model_scope = "linear",
      sim_col = "PeakFrequency",
      target_term = "contrast_num_c",
      true_beta = beta_freq,
      effect_source = "pilot_and_manifest_conservative",
      linear_nuisance_beta = 0,
      default_outcome_mean = 58.1,
      subject_dropout_rate = 0,
      output_prefix = "GCP_power_analysis_frequency",
      plot_title = "Power Analysis: Gamma Peak Frequency (subject-level condition means)"
    ),
    partition = partition_freq,
    scenario_df = build_power_scenarios(partition_freq, "residual"),
    nsim = nsim_subject,
    subject_breaks = subject_breaks,
    parallel_workers = parallel_workers,
    plot_only = plot_only
  )
  message("H5 | pessimistic scenario")
  pes_freq <- run_subject_level_power(
    config = list(
      label = "Gamma peak frequency (linear contrast)",
      hypothesis = "H5",
      pilot_outcome = "PeakFrequency",
      model_scope = "linear",
      sim_col = "PeakFrequency",
      target_term = "contrast_num_c",
      true_beta = beta_freq,
      effect_source = "pilot_and_manifest_conservative",
      linear_nuisance_beta = 0,
      default_outcome_mean = 58.1,
      output_prefix = "GCP_power_analysis_frequency_pessimistic",
      plot_title = "Power Analysis: Gamma Peak Frequency [pessimistic]"
    ),
    partition = partition_freq,
    scenario_df = build_pessimistic_scenario(partition_freq),
    nsim = nsim_subject,
    subject_breaks = subject_breaks,
    parallel_workers = parallel_workers,
    plot_only = plot_only,
    save_figures = FALSE
  )

  message("=== H6: Quadratic gamma peak power (contrast_num_c2) ===")
  message("H6 | residual sensitivity (3 residual-SD scenarios)")
  res_power <- run_subject_level_power(
    config = list(
      label = "Gamma peak power (quadratic contrast)",
      hypothesis = "H6",
      pilot_outcome = "PeakAmplitude",
      model_scope = "quadratic",
      sim_col = "PeakAmplitude",
      target_term = "contrast_num_c2",
      true_beta = beta_quad,
      effect_source = "pilot_and_manifest_conservative",
      linear_nuisance_beta = beta_linear_nuisance,
      default_outcome_mean = 3.29,
      subject_dropout_rate = 0,
      output_prefix = "GCP_power_analysis_power",
      plot_title = "Power Analysis: Gamma Peak Power (subject-level condition means)"
    ),
    partition = partition_power,
    scenario_df = build_power_scenarios(partition_power, "residual"),
    nsim = nsim_subject,
    subject_breaks = subject_breaks,
    parallel_workers = parallel_workers,
    plot_only = plot_only
  )
  message("H6 | pessimistic scenario")
  pes_power <- run_subject_level_power(
    config = list(
      label = "Gamma peak power (quadratic contrast)",
      hypothesis = "H6",
      pilot_outcome = "PeakAmplitude",
      model_scope = "quadratic",
      sim_col = "PeakAmplitude",
      target_term = "contrast_num_c2",
      true_beta = beta_quad,
      effect_source = "pilot_and_manifest_conservative",
      linear_nuisance_beta = beta_linear_nuisance,
      default_outcome_mean = 3.29,
      output_prefix = "GCP_power_analysis_power_pessimistic",
      plot_title = "Power Analysis: Gamma Peak Power [pessimistic]"
    ),
    partition = partition_power,
    scenario_df = build_pessimistic_scenario(partition_power),
    nsim = nsim_subject,
    subject_breaks = subject_breaks,
    parallel_workers = parallel_workers,
    plot_only = plot_only,
    save_figures = FALSE
  )

  message("=== H7: Gamma power x microsaccade interaction (trial-level, simr) ===")
  interaction_res <- run_interaction_power(
    gcp_root = gcp_root,
    nsim = nsim_interaction,
    subject_breaks = subject_breaks,
    power_target = power_target,
    plot_only = plot_only
  )

  summary_rows <- rbind(
    cbind(res_freq$summary_row, scenario = "residual_sensitivity"),
    data.frame(
      hypothesis = "H5",
      model = "Gamma peak frequency (linear contrast)",
      N_min_for_90 = n_for_target_power(pes_freq$power_df, power_target),
      power_at_max_N = pes_freq$power_df$power[which.max(pes_freq$power_df$n_subjects)],
      scenario = "pessimistic",
      stringsAsFactors = FALSE
    ),
    cbind(res_power$summary_row, scenario = "residual_sensitivity"),
    data.frame(
      hypothesis = "H6",
      model = "Gamma peak power (quadratic contrast)",
      N_min_for_90 = n_for_target_power(pes_power$power_df, power_target),
      power_at_max_N = pes_power$power_df$power[which.max(pes_power$power_df$n_subjects)],
      scenario = "pessimistic",
      stringsAsFactors = FALSE
    ),
    cbind(interaction_res$summary_sensitivity, scenario = "residual_sensitivity"),
    cbind(interaction_res$summary_pessimistic, scenario = "pessimistic")
  )

  pessimistic_ns <- summary_rows$N_min_for_90[summary_rows$scenario == "pessimistic"]
  binding <- summary_rows[summary_rows$scenario == "pessimistic", , drop = FALSE]
  binding <- binding[is.na(binding$N_min_for_90) | binding$N_min_for_90 == max(subject_breaks, na.rm = TRUE), , drop = FALSE]
  if (nrow(binding) > 0) {
    message(
      "Note: target power not reached within subject_breaks for: ",
      paste(unique(binding$hypothesis), collapse = ", ")
    )
  }
  recommended_N <- if (all(is.na(pessimistic_ns))) {
    max(subject_breaks)
  } else {
    max(pessimistic_ns, na.rm = TRUE)
  }

  pilot_ged_exclusion_rate <- 0.30
  recommendation <- data.frame(
    recommended_N_analyzable = recommended_N,
    power_target = power_target,
    rule = "max N_min_for_90 across H5, H6, H7 under pessimistic variance scenarios",
    recruitment_rule = "recruit until recommended_N_analyzable with valid GED (no simulated dropout)",
    subject_dropout_in_simulation = 0,
    pilot_ged_exclusion_rate = pilot_ged_exclusion_rate,
    pilot_ged_exclusion_n = 3L,
    pilot_ged_screened_n = 10L,
    enrollment_if_30pct_exclusion = as.integer(ceiling(recommended_N / (1 - pilot_ged_exclusion_rate))),
    beta_frequency = beta_freq,
    beta_quadratic_power = beta_quad,
    beta_interaction = interaction_res$beta_interaction,
    effect_note = "Betas from pilot/manifest conservative rule; Table 1 Cohen's d not used in simulations",
    pilot_n_frequency = partition_freq$pilot_n_subjects,
    pilot_n_power = partition_power$pilot_n_subjects,
    stringsAsFactors = FALSE
  )

  utils::write.csv(summary_rows, file.path(data_dir, "GCP_power_analysis_summary.csv"), row.names = FALSE)
  utils::write.csv(recommendation, file.path(data_dir, "GCP_power_analysis_recommended_N.csv"), row.names = FALSE)
  utils::write.csv(pes_freq$power_df, file.path(data_dir, "GCP_power_analysis_frequency_pessimistic_curve.csv"), row.names = FALSE)
  utils::write.csv(pes_power$power_df, file.path(data_dir, "GCP_power_analysis_power_pessimistic_curve.csv"), row.names = FALSE)
  message("Wrote summary CSVs to: ", data_dir)

  message(sprintf("\n=== Recommended sample size: N = %d (target power = %.0f%%) ===", recommended_N, 100 * power_target))
  print(recommendation)
  print(summary_rows)

  invisible(list(
    recommendation = recommendation,
    summary = summary_rows,
    frequency = res_freq,
    power = res_power,
    interaction = interaction_res
  ))
}

extract_powercurve_df <- function(pc, n_subjects_vec) {
  do.call(rbind, lapply(seq_along(pc$ps), function(i) {
    z <- pc$ps[[i]]
    ci <- suppressWarnings(stats::confint(z))
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

run_interaction_power <- function(
    gcp_root,
    nsim = 500L,
    subject_breaks = seq(20, 100, by = 10),
    power_target = 0.90,
    plot_only = FALSE) {
  pilot_path <- file.path(gcp_root, "data", "pilot_stats", "pilot_interaction_power_parameters.csv")
  if (!file.exists(pilot_path)) {
    stop("Interaction pilot parameters not found: ", pilot_path)
  }
  pilot <- utils::read.csv(pilot_path, stringsAsFactors = FALSE)
  if (nrow(pilot) != 1L) {
    stop("Expected one row in pilot_interaction_power_parameters.csv")
  }

  contrast_levels <- c("25", "50", "75", "100")
  trials_per_condition <- 160L
  seed <- 123L

  true_beta_interaction <- as.numeric(pilot$beta_interaction)
  true_beta_contrast <- as.numeric(pilot$beta_contrast)
  true_beta_microsaccade <- as.numeric(pilot$beta_microsaccade_scaled)
  outcome_mean <- as.numeric(pilot$outcome_mean)
  microsaccade_mean <- as.numeric(pilot$microsaccade_grand_mean)
  baseline_random_intercept_sd <- as.numeric(pilot$random_intercept_sd)
  baseline_residual_levels <- c(
    as.numeric(pilot$residual_sd_q025),
    as.numeric(pilot$residual_sd_q50),
    as.numeric(pilot$residual_sd_q975)
  )
  pessimistic_residual <- as.numeric(pilot$residual_sd_q975)
  ms_subject_sd <- as.numeric(pilot$ms_random_intercept_sd)
  ms_residual_sd <- as.numeric(pilot$ms_residual_sd)
  ms_effect_contrast <- as.numeric(pilot$ms_beta_contrast_raw)

  make_template <- function(n_subjects) {
    Subject <- factor(rep(seq_len(n_subjects), each = length(contrast_levels) * trials_per_condition))
    contrast <- factor(
      rep(rep(contrast_levels, each = trials_per_condition), times = n_subjects),
      levels = contrast_levels, ordered = TRUE
    )
    contrast_num <- as.numeric(as.character(contrast))
    data.frame(
      Subject = Subject,
      contrast_num_c = as.numeric(scale(contrast_num, center = TRUE, scale = TRUE))
    )
  }

  simulate_microsaccade <- function(dat) {
    n_sub <- nlevels(dat$Subject)
    ms_intercept_subj <- stats::rnorm(n_sub, mean = 0, sd = ms_subject_sd)
    ms_raw <- microsaccade_mean +
      ms_effect_contrast * dat$contrast_num_c +
      ms_intercept_subj[dat$Subject] +
      stats::rnorm(nrow(dat), mean = 0, sd = ms_residual_sd)
    as.numeric(scale(ms_raw, center = TRUE, scale = TRUE))
  }

  simulate_response <- function(dat, residual_sd) {
    n_sub <- nlevels(dat$Subject)
    u0 <- stats::rnorm(n_sub, mean = 0, sd = baseline_random_intercept_sd)
    x <- dat$contrast_num_c
    m <- dat$microsaccade_c
    y_mu <- outcome_mean +
      u0[dat$Subject] +
      true_beta_contrast * x +
      true_beta_microsaccade * m +
      true_beta_interaction * x * m
    y_mu + stats::rnorm(nrow(dat), mean = 0, sd = residual_sd)
  }

  run_simr_scenarios <- function(residual_levels, labels, output_prefix, plot_title) {
    message(sprintf("  Interaction block: %s (%d scenario(s), nsim=%d)", output_prefix, length(labels), nsim))
    set.seed(seed)
    template_dat <- make_template(max(subject_breaks))
    template_dat$microsaccade_c <- simulate_microsaccade(template_dat)
    template_dat$gamma_power <- simulate_response(template_dat, baseline_residual_levels[2])

    message("  Fitting template mixed model for simr...")
    base_model <- lmer(
      gamma_power ~ contrast_num_c * microsaccade_c + (1 | Subject),
      data = template_dat,
      REML = FALSE,
      control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
    )
    lme4::fixef(base_model)["(Intercept)"] <- outcome_mean
    lme4::fixef(base_model)["contrast_num_c"] <- true_beta_contrast
    lme4::fixef(base_model)["microsaccade_c"] <- true_beta_microsaccade
    lme4::fixef(base_model)["contrast_num_c:microsaccade_c"] <- true_beta_interaction
    lme4::VarCorr(base_model) <- lme4::VarCorr(
      lmer(gamma_power ~ contrast_num_c * microsaccade_c + (1 | Subject), data = template_dat, REML = FALSE)
    )

    scenario_df <- data.frame(
      scenario_label = labels,
      scenario_display = sprintf("Residual %.3f", residual_levels),
      residual_sd = residual_levels,
      stringsAsFactors = FALSE
    )

    power_rows <- lapply(seq_len(nrow(scenario_df)), function(i) {
      sc <- scenario_df[i, , drop = FALSE]
      message(sprintf("  Scenario %s | %s | residual SD = %.3f", sc$scenario_label, sc$scenario_display, sc$residual_sd))
      model_sc <- base_model
      stats::sigma(model_sc) <- sc$residual_sd
      pc <- suppressWarnings(powerCurve(
        model_sc,
        test = fixed("contrast_num_c:microsaccade_c"),
        along = "Subject",
        breaks = subject_breaks,
        nsim = nsim,
        progress = TRUE
      ))
      out <- extract_powercurve_df(pc, subject_breaks)
      for (j in seq_len(nrow(out))) {
        message(sprintf(
          "    N=%3d | power=%.3f | nsim=%d",
          out$n_subjects[j], out$power[j], out$nsim[j]
        ))
      }
      out$scenario_label <- sc$scenario_label
      out$scenario_display <- sc$scenario_display
      out$residual_sd <- sc$residual_sd
      out
    })
    power_df <- do.call(rbind, power_rows)

    figure_dir <- file.path(gcp_root, "figures", "power_analysis")
    data_dir <- file.path(gcp_root, "data", "power_analysis")
    dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
    curve_path <- file.path(data_dir, paste0(output_prefix, "_curve.csv"))
    utils::write.csv(power_df, curve_path, row.names = FALSE)
    message(sprintf("  Wrote interaction curve: %s", curve_path))

    if (nrow(scenario_df) > 1L) {
      message(sprintf("  Saving interaction figures: %s", output_prefix))
      save_power_figures(
        power_df = power_df,
        scenario_levels = scenario_df$scenario_display,
        subject_breaks = subject_breaks,
        strict_power_target = power_target,
        figure_output_dir = figure_dir,
        output_prefix = output_prefix,
        plot_title = plot_title,
        y_legend_title = "Residuals SD"
      )
    }

    list(power_df = power_df, scenario_df = scenario_df)
  }

  if (!plot_only) {
    message("H7 | residual sensitivity (3 residual-SD scenarios)")
    sens_labels <- c("res_low", "res_median", "res_high")
    sens <- run_simr_scenarios(
      baseline_residual_levels, sens_labels,
      "GCP_power_analysis_interaction",
      "Power Analysis: Gamma Power x Microsaccade Interaction (trial-level)"
    )
    message("H7 | pessimistic scenario")
    pes <- run_simr_scenarios(
      c(pessimistic_residual), c("pessimistic"),
      "GCP_power_analysis_interaction_pessimistic",
      "Power Analysis: Gamma Power x Microsaccade Interaction [pessimistic]"
    )
    power_df <- sens$power_df
    pessimistic_df <- pes$power_df
  } else {
    power_df <- utils::read.csv(
      file.path(gcp_root, "data", "power_analysis", "GCP_power_analysis_interaction_curve.csv"),
      stringsAsFactors = FALSE
    )
    pessimistic_df <- utils::read.csv(
      file.path(gcp_root, "data", "power_analysis", "GCP_power_analysis_interaction_pessimistic_curve.csv"),
      stringsAsFactors = FALSE
    )
  }

  high_label <- sprintf("Residual %.3f", baseline_residual_levels[3])
  summary_sensitivity <- data.frame(
    hypothesis = "H7",
    model = "Gamma power x microsaccade interaction",
    N_min_for_90 = n_for_target_power(
      power_df[power_df$scenario_display == high_label, , drop = FALSE], power_target
    ),
    power_at_max_N = {
      sub <- power_df[power_df$scenario_display == high_label, , drop = FALSE]
      sub$power[which.max(sub$n_subjects)]
    },
    stringsAsFactors = FALSE
  )
  summary_pessimistic <- data.frame(
    hypothesis = "H7",
    model = "Gamma power x microsaccade interaction",
    N_min_for_90 = n_for_target_power(pessimistic_df, power_target),
    power_at_max_N = pessimistic_df$power[which.max(pessimistic_df$n_subjects)],
    stringsAsFactors = FALSE
  )

  list(
    power_df = power_df,
    pessimistic_df = pessimistic_df,
    summary_sensitivity = summary_sensitivity,
    summary_pessimistic = summary_pessimistic,
    beta_interaction = true_beta_interaction
  )
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  plot_only <- "--plot-only" %in% args
  nsim_subject <- 1000L
  nsim_interaction <- 500L
  hit <- grep("^--nsim=", args, value = TRUE)
  if (length(hit) > 0) {
    nsim_subject <- as.integer(sub("^--nsim=", "", hit[1]))
    nsim_interaction <- max(200L, nsim_subject %/% 2L)
  }
  run_power_analysis(
    nsim_subject = nsim_subject,
    nsim_interaction = nsim_interaction,
    plot_only = plot_only
  )
}
