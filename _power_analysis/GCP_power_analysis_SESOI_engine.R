## Shared SESOI subject-level power simulation engine.

source_power_helpers <- function() {
  helper_path <- file.path(getwd(), "_power_analysis", "GCP_power_analysis_helpers.R")
  if (!file.exists(helper_path)) {
    cmd <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", cmd, value = TRUE)
    if (length(file_arg) > 0) {
      helper_path <- file.path(dirname(sub("^--file=", "", file_arg[1])), "GCP_power_analysis_helpers.R")
    }
  }
  if (!file.exists(helper_path)) {
    stop("Could not find GCP_power_analysis_helpers.R")
  }
  source(helper_path)
  helper_path
}

power_binomial_se <- function(successes, n) {
  if (n <= 0) {
    return(NA_real_)
  }
  p <- successes / n
  sqrt(p * (1 - p) / n)
}

run_sesoi_power <- function(
    config,
    plot_only = FALSE,
    nsim = 1000,
    subject_breaks = seq(20, 100, by = 10),
    parallel_workers = 8,
    sensitivity_axes = c("residual", "random_effects", "fixed_multiplier")) {
  helper_path <- source_power_helpers()
  suppressPackageStartupMessages(library(lme4))
  if (requireNamespace("lmerTest", quietly = TRUE)) {
    suppressPackageStartupMessages(library(lmerTest))
  }

  gcp_root <- resolve_gcp_root()
  partition <- load_variance_partition_summary(gcp_root, config$pilot_outcome, config$model_scope)
  formula_list <- build_subject_level_formula_list(config$sim_col, model_scope = config$model_scope)

  manifest <- load_manifest_effects(
    gcp_root = gcp_root,
    pilot_outcome = config$pilot_outcome,
    target_term = config$target_term,
    true_beta_multiplier = config$true_beta_multiplier %||% 1.5,
    fallback_sesoi_beta = config$sesoi_beta,
    fallback_true_beta = config$true_beta
  )
  config$sesoi_beta <- manifest$sesoi_beta
  config$true_beta <- manifest$true_beta
  config$manifest_source <- manifest$source

  outcome_mean <- partition$point_estimates$outcome_mean
  if (is.null(outcome_mean) || !is.finite(outcome_mean)) {
    outcome_mean <- config$default_outcome_mean
  }
  linear_nuisance_beta <- partition$point_estimates$linear_nuisance_beta
  if (is.null(linear_nuisance_beta) || !is.finite(linear_nuisance_beta)) {
    linear_nuisance_beta <- config$linear_nuisance_beta %||% 0
  }
  config$linear_nuisance_beta <- linear_nuisance_beta

  message(sprintf(
    "Manifest SESOI beta=%.4f | true beta=%.4f | source=%s",
    config$sesoi_beta, config$true_beta, config$manifest_source
  ))

  results <- list()
  for (axis in sensitivity_axes) {
    message(sprintf("\n=== %s | sensitivity axis: %s ===", config$label, axis))
    results[[axis]] <- run_sesoi_power_axis(
      config = config,
      partition = partition,
      formula_list = formula_list,
      outcome_mean = outcome_mean,
      sensitivity_axis = axis,
      plot_only = plot_only,
      nsim = nsim,
      subject_breaks = subject_breaks,
      parallel_workers = parallel_workers,
      helper_path = helper_path,
      gcp_root = gcp_root
    )
  }
  results
}

run_sesoi_power_axis <- function(
    config,
    partition,
    formula_list,
    outcome_mean,
    sensitivity_axis,
    plot_only,
    nsim,
    subject_breaks,
    parallel_workers,
    helper_path,
    gcp_root) {
  seed <- 123
  sesoi_decision_alpha <- 0.05
  detection_alpha <- 0.05
  strict_power_target <- 0.90
  contrast_levels <- c("25", "50", "75", "100")
  subject_dropout_rate <- config$subject_dropout_rate %||% 0.10
  parallel_round_chunk_nsim <- 10L

  scenario_df <- build_power_scenarios(partition, sensitivity_axis)
  output_prefix <- paste0(config$output_prefix_base, "_", sensitivity_axis)
  plot_title <- paste0(config$plot_title_base, " [", sensitivity_axis, "]")
  figure_output_dir <- file.path(gcp_root, "figures", "power_analysis")
  data_output_dir <- file.path(gcp_root, "data", "power_analysis")
  curve_csv_path <- file.path(data_output_dir, paste0(output_prefix, "_curve.csv"))

  dir.create(figure_output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(data_output_dir, recursive = TRUE, showWarnings = FALSE)

  message(
    "Pilot N=", partition$pilot_n_subjects,
    " | analysis unit: one condition mean per subject (",
    SUBJECT_CONDITIONS_PER_SUBJECT, " rows/subject)"
  )
  message(
    "Identified RE structure: RI=",
    partition$includes_random_intercept,
    " RS=", partition$includes_random_slope,
    if (config$model_scope == "quadratic") paste0(" RQS=", partition$includes_random_quadratic_slope) else "",
    " | formula: ", partition$lmer_formula
  )
  print(scenario_df[, c("scenario_label", "random_intercept_sd", "random_slope_sd",
                        "random_quadratic_slope_sd", "residual_sd")], row.names = FALSE)

  estimate_power_chunk <- function(
      chunk_nsim,
      n_subjects,
      formula_list,
      sesoi_decision_alpha,
      detection_alpha,
      beta_raw,
      true_beta,
      outcome_mean,
      random_intercept_sd,
      random_slope_sd,
      random_quadratic_slope_sd,
      residual_sd,
      linear_nuisance_beta,
      use_random_slope,
      use_random_quadratic_slope,
      subject_dropout_rate,
      sim_col,
      target_term,
      model_scope) {
    valid_fits <- 0L
    sesoi_successes <- 0L
    detection_successes <- 0L
    sign_errors <- 0L
    type_m_sum <- 0
    type_m_count <- 0L
    singular_count <- 0L
    convergence_warn_count <- 0L
    rs_sd <- if (use_random_slope) random_slope_sd else 0
    rqs_sd <- if (use_random_quadratic_slope) random_quadratic_slope_sd else 0

    for (i in seq_len(chunk_nsim)) {
      dat <- make_subject_level_design(n_subjects, contrast_levels)
      n_subject_levels <- nlevels(dat$Subject)
      random_intercepts <- stats::rnorm(n_subject_levels, mean = 0, sd = random_intercept_sd)
      random_slopes <- stats::rnorm(n_subject_levels, mean = 0, sd = rs_sd)
      random_quadratic_slopes <- stats::rnorm(n_subject_levels, mean = 0, sd = rqs_sd)
      x <- dat$contrast_num_c
      x2 <- dat$contrast_num_c2
      mu <- outcome_mean + random_intercepts[dat$Subject]
      if (model_scope == "quadratic") {
        mu <- mu +
          random_slopes[dat$Subject] * x +
          random_quadratic_slopes[dat$Subject] * x2 +
          linear_nuisance_beta * x +
          true_beta * x2
      } else {
        mu <- mu + random_slopes[dat$Subject] * x + true_beta * x
      }
      dat[[sim_col]] <- mu + stats::rnorm(nrow(dat), mean = 0, sd = residual_sd)

      dat <- apply_subject_dropout(dat, subject_dropout_rate)
      if (!subject_level_fit_eligible(dat)) {
        next
      }

      fit_info <- fit_lmer_with_fallbacks(formula_list, dat)
      if (is.null(fit_info)) {
        next
      }
      fit_full <- fit_info$fit
      infer_out <- inference_from_fit(
        fit = fit_full,
        term = target_term,
        beta_raw = beta_raw,
        sesoi_decision_alpha = sesoi_decision_alpha,
        true_beta = true_beta,
        detection_alpha = detection_alpha
      )
      if (!infer_out$valid) {
        next
      }

      valid_fits <- valid_fits + 1L
      singular_count <- singular_count + as.integer(fit_info$singular)
      conv_msgs <- fit_full@optinfo$conv$lme4$messages
      convergence_warn_count <- convergence_warn_count + as.integer(length(conv_msgs) > 0)
      sesoi_successes <- sesoi_successes + as.integer(infer_out$sesoi_success)
      detection_successes <- detection_successes + as.integer(infer_out$detection_success)
      sign_errors <- sign_errors + as.integer(infer_out$sign_error)
      if (isTRUE(infer_out$sesoi_success)) {
        type_m_sum <- type_m_sum + abs(infer_out$est / true_beta)
        type_m_count <- type_m_count + 1L
      }
    }

    list(
      chunk_nsim = as.integer(chunk_nsim),
      valid_fits = as.integer(valid_fits),
      sesoi_successes = as.integer(sesoi_successes),
      detection_successes = as.integer(detection_successes),
      sign_errors = as.integer(sign_errors),
      type_m_sum = as.numeric(type_m_sum),
      type_m_count = as.integer(type_m_count),
      singular_count = as.integer(singular_count),
      convergence_warn_count = as.integer(convergence_warn_count)
    )
  }

  if (!plot_only) {
    set.seed(seed)
    cl <- parallel::makeCluster(as.integer(parallel_workers), type = "PSOCK")
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterEvalQ(cl, {
      suppressPackageStartupMessages(library(lme4))
      if (requireNamespace("lmerTest", quietly = TRUE)) {
        suppressPackageStartupMessages(library(lmerTest))
      }
    })
    parallel::clusterExport(cl, "helper_path", envir = environment())
    parallel::clusterEvalQ(cl, source(helper_path))
    parallel::clusterSetRNGStream(cl, iseed = seed)
    parallel::clusterExport(
      cl,
      varlist = c(
        "estimate_power_chunk", "contrast_levels", "make_subject_level_design",
        "apply_subject_dropout", "fit_lmer_with_fallbacks", "inference_from_fit",
        "subject_level_fit_eligible", "SUBJECT_CONDITIONS_PER_SUBJECT"
      ),
      envir = environment()
    )

    scenario_order <- scenario_df$scenario_label
    scenario_rows <- lapply(seq_len(nrow(scenario_df)), function(idx) {
      scenario <- scenario_df[idx, , drop = FALSE]
      n_rows <- lapply(subject_breaks, function(n_subjects) {
        cat(sprintf(
          "[%s] START axis=%s scenario=%s N=%d",
          format(Sys.time(), "%H:%M:%S"), sensitivity_axis, scenario$scenario_label, n_subjects
        ))
        flush.console()

        total_nsim <- 0L
        total_valid_fits <- 0L
        total_sesoi_successes <- 0L
        total_detection_successes <- 0L
        total_sign_errors <- 0L
        total_type_m_sum <- 0
        total_type_m_count <- 0L
        total_singular <- 0L
        total_convergence_warn <- 0L

        while (total_nsim < nsim) {
          remaining <- nsim - total_nsim
          workers_this_round <- min(length(cl), ceiling(remaining / parallel_round_chunk_nsim))
          workers_this_round <- max(workers_this_round, 1L)
          chunk_sizes <- rep(parallel_round_chunk_nsim, workers_this_round)
          total_chunk <- sum(chunk_sizes)
          if (total_chunk > remaining) {
            overflow <- total_chunk - remaining
            for (k in seq_len(overflow)) {
              idx2 <- ((k - 1) %% workers_this_round) + 1L
              chunk_sizes[idx2] <- chunk_sizes[idx2] - 1L
            }
          }
          chunk_sizes <- chunk_sizes[chunk_sizes > 0L]

          chunk_results <- parallel::parLapply(
            cl = cl,
            X = chunk_sizes,
            fun = estimate_power_chunk,
            n_subjects = n_subjects,
            formula_list = formula_list,
            sesoi_decision_alpha = sesoi_decision_alpha,
            detection_alpha = detection_alpha,
            beta_raw = config$sesoi_beta,
            true_beta = config$true_beta,
            outcome_mean = outcome_mean,
            random_intercept_sd = scenario$random_intercept_sd,
            random_slope_sd = scenario$random_slope_sd,
            random_quadratic_slope_sd = scenario$random_quadratic_slope_sd,
            residual_sd = scenario$residual_sd,
            linear_nuisance_beta = config$linear_nuisance_beta %||% 0,
            use_random_slope = partition$includes_random_slope,
            use_random_quadratic_slope = partition$includes_random_quadratic_slope,
            subject_dropout_rate = subject_dropout_rate,
            sim_col = config$sim_col,
            target_term = config$target_term,
            model_scope = config$model_scope
          )

          total_nsim <- total_nsim + sum(vapply(chunk_results, function(x) x$chunk_nsim, integer(1)))
          total_valid_fits <- total_valid_fits + sum(vapply(chunk_results, function(x) x$valid_fits, integer(1)))
          total_sesoi_successes <- total_sesoi_successes + sum(vapply(chunk_results, function(x) x$sesoi_successes, integer(1)))
          total_detection_successes <- total_detection_successes + sum(vapply(chunk_results, function(x) x$detection_successes, integer(1)))
          total_sign_errors <- total_sign_errors + sum(vapply(chunk_results, function(x) x$sign_errors, integer(1)))
          total_type_m_sum <- total_type_m_sum + sum(vapply(chunk_results, function(x) x$type_m_sum, numeric(1)))
          total_type_m_count <- total_type_m_count + sum(vapply(chunk_results, function(x) x$type_m_count, integer(1)))
          total_singular <- total_singular + sum(vapply(chunk_results, function(x) x$singular_count, integer(1)))
          total_convergence_warn <- total_convergence_warn + sum(vapply(chunk_results, function(x) x$convergence_warn_count, integer(1)))
        }

        power_sesoi_cond <- if (total_valid_fits > 0) total_sesoi_successes / total_valid_fits else NA_real_
        power_sesoi_uncond <- if (total_nsim > 0) total_sesoi_successes / total_nsim else NA_real_
        power_detection_cond <- if (total_valid_fits > 0) total_detection_successes / total_valid_fits else NA_real_
        power_detection_uncond <- if (total_nsim > 0) total_detection_successes / total_nsim else NA_real_
        cat(sprintf(
          " | SESOI uncond=%.3f detect uncond=%.3f | fit=%.0f%%\n",
          power_sesoi_uncond, power_detection_uncond, 100 * total_valid_fits / total_nsim
        ))
        flush.console()

        se_sesoi_uncond <- power_binomial_se(total_sesoi_successes, total_nsim)
        se_detection_uncond <- power_binomial_se(total_detection_successes, total_nsim)

        data.frame(
          scenario_label = scenario$scenario_label,
          scenario_display = scenario$scenario_display,
          sensitivity_axis = sensitivity_axis,
          varied_component = scenario$varied_component,
          random_intercept_sd = scenario$random_intercept_sd,
          random_slope_sd = scenario$random_slope_sd,
          random_quadratic_slope_sd = scenario$random_quadratic_slope_sd,
          residual_sd = scenario$residual_sd,
          n_subjects = n_subjects,
          sesoi_beta = config$sesoi_beta,
          true_beta = config$true_beta,
          power_sesoi_conditional = power_sesoi_cond,
          power_sesoi_unconditional = power_sesoi_uncond,
          power_detection_conditional = power_detection_cond,
          power_detection_unconditional = power_detection_uncond,
          lower_sesoi_unconditional = if (is.finite(se_sesoi_uncond)) {
            pmax(0, power_sesoi_uncond - 1.96 * se_sesoi_uncond)
          } else {
            NA_real_
          },
          upper_sesoi_unconditional = if (is.finite(se_sesoi_uncond)) {
            pmin(1, power_sesoi_uncond + 1.96 * se_sesoi_uncond)
          } else {
            NA_real_
          },
          lower_detection_unconditional = if (is.finite(se_detection_uncond)) {
            pmax(0, power_detection_uncond - 1.96 * se_detection_uncond)
          } else {
            NA_real_
          },
          upper_detection_unconditional = if (is.finite(se_detection_uncond)) {
            pmin(1, power_detection_uncond + 1.96 * se_detection_uncond)
          } else {
            NA_real_
          },
          power = power_sesoi_uncond,
          lower = if (is.finite(se_sesoi_uncond)) pmax(0, power_sesoi_uncond - 1.96 * se_sesoi_uncond) else NA_real_,
          upper = if (is.finite(se_sesoi_uncond)) pmin(1, power_sesoi_uncond + 1.96 * se_sesoi_uncond) else NA_real_,
          nsim = total_nsim,
          valid_fits = total_valid_fits,
          fit_success_rate = if (total_nsim > 0) total_valid_fits / total_nsim else NA_real_,
          type_s = if (total_valid_fits > 0) total_sign_errors / total_valid_fits else NA_real_,
          type_m = if (total_type_m_count > 0) total_type_m_sum / total_type_m_count else NA_real_,
          singular_rate = if (total_valid_fits > 0) total_singular / total_valid_fits else NA_real_,
          convergence_warn_rate = if (total_valid_fits > 0) total_convergence_warn / total_valid_fits else NA_real_,
          analysis_unit = "subject_level_condition_mean",
          observations_per_subject = SUBJECT_CONDITIONS_PER_SUBJECT,
          pilot_n_subjects = partition$pilot_n_subjects,
          model_scope = config$model_scope,
          lmer_formula_pilot = partition$lmer_formula,
          manifest_source = config$manifest_source,
          stringsAsFactors = FALSE
        )
      })
      do.call(rbind, n_rows)
    })

    power_df <- do.call(rbind, scenario_rows)
    power_df$scenario_label <- factor(power_df$scenario_label, levels = scenario_order, ordered = TRUE)
    power_df$meets_target_90 <- power_df$power_sesoi_unconditional >= strict_power_target
    utils::write.csv(power_df, curve_csv_path, row.names = FALSE)
  } else {
    if (!file.exists(curve_csv_path)) {
      stop("plot_only=TRUE requested but curve CSV not found: ", curve_csv_path)
    }
    power_df <- utils::read.csv(curve_csv_path, stringsAsFactors = FALSE)
    if (!"meets_target_90" %in% names(power_df)) {
      primary <- if ("power_sesoi_unconditional" %in% names(power_df)) {
        power_df$power_sesoi_unconditional
      } else {
        power_df$power
      }
      power_df$meets_target_90 <- primary >= strict_power_target
    }
  }

  summary_rows <- do.call(rbind, lapply(split(power_df, power_df$scenario_label), function(df) {
    df <- df[order(df$n_subjects), , drop = FALSE]
    hit <- df[df$meets_target_90, "n_subjects", drop = TRUE]
    primary_col <- if ("power_sesoi_unconditional" %in% names(df)) {
      "power_sesoi_unconditional"
    } else {
      "power"
    }
    detect_col <- if ("power_detection_unconditional" %in% names(df)) {
      "power_detection_unconditional"
    } else {
      primary_col
    }
    data.frame(
      scenario_label = as.character(df$scenario_label[1]),
      sensitivity_axis = sensitivity_axis,
      varied_component = as.character(df$varied_component[1]),
      N_min_for_90_sesoi = if (length(hit) > 0) min(hit) else NA_integer_,
      sesoi_uncond_at_max_N = df[[primary_col]][which.max(df$n_subjects)],
      detection_uncond_at_max_N = df[[detect_col]][which.max(df$n_subjects)],
      fit_success_rate_at_max_N = df$fit_success_rate[which.max(df$n_subjects)],
      type_s_at_max_N = df$type_s[which.max(df$n_subjects)],
      type_m_at_max_N = df$type_m[which.max(df$n_subjects)],
      singular_rate_at_max_N = df$singular_rate[which.max(df$n_subjects)],
      convergence_warn_rate_at_max_N = df$convergence_warn_rate[which.max(df$n_subjects)],
      monotonic_non_decreasing = all(diff(df[[primary_col]]) >= -0.02),
      stringsAsFactors = FALSE
    )
  }))
  utils::write.csv(
    summary_rows,
    file.path(data_output_dir, paste0(output_prefix, "_summary.csv")),
    row.names = FALSE
  )

  save_power_figures(
    power_df = power_df,
    scenario_levels = unique(scenario_df$scenario_display),
    subject_breaks = subject_breaks,
    strict_power_target = strict_power_target,
    figure_output_dir = figure_output_dir,
    output_prefix = output_prefix,
    plot_title = plot_title,
    sensitivity_axis = sensitivity_axis
  )

  power_df
}
