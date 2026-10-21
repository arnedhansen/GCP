## Monte Carlo detection-power engine for subject-level mixed models.

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

run_subject_level_power <- function(
    config,
    partition,
    scenario_df,
    nsim = 1000,
    subject_breaks = seq(20, 100, by = 10),
    parallel_workers = 8,
    plot_only = FALSE,
    save_figures = TRUE) {
  helper_path <- source_power_helpers()
  suppressPackageStartupMessages(library(lme4))
  if (requireNamespace("lmerTest", quietly = TRUE)) {
    suppressPackageStartupMessages(library(lmerTest))
  }

  gcp_root <- resolve_gcp_root()
  formula_list <- build_subject_level_formula_list(config$sim_col, model_scope = config$model_scope)
  outcome_mean <- partition$point_estimates$outcome_mean
  if (is.null(outcome_mean) || !is.finite(outcome_mean)) {
    outcome_mean <- config$default_outcome_mean
  }
  linear_nuisance_beta <- partition$point_estimates$linear_nuisance_beta
  if (is.null(linear_nuisance_beta) || !is.finite(linear_nuisance_beta)) {
    linear_nuisance_beta <- config$linear_nuisance_beta %||% 0
  }

  seed <- 123
  detection_alpha <- 0.05
  strict_power_target <- 0.90
  contrast_levels <- c("25", "50", "75", "100")
  parallel_round_chunk_nsim <- 10L

  figure_output_dir <- file.path(gcp_root, "figures", "power_analysis")
  data_output_dir <- file.path(gcp_root, "data", "power_analysis")
  curve_csv_path <- file.path(data_output_dir, paste0(config$output_prefix, "_curve.csv"))
  dir.create(figure_output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(data_output_dir, recursive = TRUE, showWarnings = FALSE)

  message(sprintf(
    "%s | true beta=%.4f | outcome mean=%.3f",
    config$label, config$true_beta, outcome_mean
  ))
  if (plot_only) {
    message(sprintf("  [plot_only] reading cached curve: %s", curve_csv_path))
  } else {
    message(sprintf(
      "  Starting: %s | %d scenario(s) x %d sample sizes | nsim=%d | workers=%d",
      config$output_prefix, nrow(scenario_df), length(subject_breaks), nsim, parallel_workers
    ))
  }

  estimate_power_chunk <- function(
      chunk_nsim, n_subjects, formula_list, detection_alpha,
      true_beta, outcome_mean, random_intercept_sd, random_slope_sd,
      random_quadratic_slope_sd, residual_sd, linear_nuisance_beta,
      use_random_slope, use_random_quadratic_slope,
      sim_col, target_term, model_scope) {
    valid_fits <- 0L
    detection_successes <- 0L
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
      if (!subject_level_fit_eligible(dat)) {
        next
      }
      fit_info <- fit_lmer_with_fallbacks(formula_list, dat)
      if (is.null(fit_info)) {
        next
      }
      det_out <- detection_from_fit(fit_info$fit, target_term, detection_alpha)
      if (!det_out$valid) {
        next
      }
      valid_fits <- valid_fits + 1L
      detection_successes <- detection_successes + as.integer(det_out$detection_success)
    }

    list(
      chunk_nsim = as.integer(chunk_nsim),
      valid_fits = as.integer(valid_fits),
      detection_successes = as.integer(detection_successes)
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
        "fit_lmer_with_fallbacks", "detection_from_fit",
        "subject_level_fit_eligible", "SUBJECT_CONDITIONS_PER_SUBJECT"
      ),
      envir = environment()
    )

    scenario_order <- scenario_df$scenario_label
    scenario_rows <- lapply(seq_len(nrow(scenario_df)), function(idx) {
      scenario <- scenario_df[idx, , drop = FALSE]
      message(sprintf(
        "  Scenario %s | %s | residual SD=%.3f",
        scenario$scenario_label, scenario$scenario_display, scenario$residual_sd
      ))
      n_rows <- lapply(subject_breaks, function(n_subjects) {
        total_nsim <- 0L
        total_valid_fits <- 0L
        total_detection_successes <- 0L
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
            cl, chunk_sizes, estimate_power_chunk,
            n_subjects = n_subjects,
            formula_list = formula_list,
            detection_alpha = detection_alpha,
            true_beta = config$true_beta,
            outcome_mean = outcome_mean,
            random_intercept_sd = scenario$random_intercept_sd,
            random_slope_sd = scenario$random_slope_sd,
            random_quadratic_slope_sd = scenario$random_quadratic_slope_sd,
            residual_sd = scenario$residual_sd,
            linear_nuisance_beta = linear_nuisance_beta,
            use_random_slope = partition$includes_random_slope,
            use_random_quadratic_slope = partition$includes_random_quadratic_slope,
            sim_col = config$sim_col,
            target_term = config$target_term,
            model_scope = config$model_scope
          )
          total_nsim <- total_nsim + sum(vapply(chunk_results, function(x) x$chunk_nsim, integer(1)))
          total_valid_fits <- total_valid_fits + sum(vapply(chunk_results, function(x) x$valid_fits, integer(1)))
          total_detection_successes <- total_detection_successes +
            sum(vapply(chunk_results, function(x) x$detection_successes, integer(1)))
        }
        power <- if (total_nsim > 0) total_detection_successes / total_nsim else NA_real_
        se <- power_binomial_se(total_detection_successes, total_nsim)
        fit_rate <- if (total_nsim > 0) total_valid_fits / total_nsim else NA_real_
        message(sprintf(
          "    N=%3d | power=%.3f | nsim=%d | valid fits=%.1f%%",
          n_subjects, power, total_nsim, 100 * fit_rate
        ))
        data.frame(
          scenario_label = scenario$scenario_label,
          scenario_display = scenario$scenario_display,
          varied_component = scenario$varied_component,
          random_intercept_sd = scenario$random_intercept_sd,
          random_slope_sd = scenario$random_slope_sd,
          random_quadratic_slope_sd = scenario$random_quadratic_slope_sd,
          residual_sd = scenario$residual_sd,
          n_subjects = n_subjects,
          true_beta = config$true_beta,
          effect_source = config$effect_source %||% "unspecified",
          power = power,
          lower = if (is.finite(se)) pmax(0, power - 1.96 * se) else NA_real_,
          upper = if (is.finite(se)) pmin(1, power + 1.96 * se) else NA_real_,
          nsim = total_nsim,
          valid_fits = total_valid_fits,
          fit_success_rate = if (total_nsim > 0) total_valid_fits / total_nsim else NA_real_,
          analysis_unit = "subject_level_condition_mean",
          model_scope = config$model_scope,
          hypothesis = config$hypothesis,
          stringsAsFactors = FALSE
        )
      })
      do.call(rbind, n_rows)
    })

    power_df <- do.call(rbind, scenario_rows)
    power_df$scenario_label <- factor(power_df$scenario_label, levels = scenario_order, ordered = TRUE)
    utils::write.csv(power_df, curve_csv_path, row.names = FALSE)
    message(sprintf("  Wrote power curve: %s", curve_csv_path))
  } else {
    if (!file.exists(curve_csv_path)) {
      stop("plot_only=TRUE but CSV not found: ", curve_csv_path)
    }
    power_df <- utils::read.csv(curve_csv_path, stringsAsFactors = FALSE)
  }

  if (nrow(scenario_df) == 1L) {
    summary_row <- data.frame(
      hypothesis = config$hypothesis,
      model = config$label,
      N_min_for_90 = n_for_target_power(power_df, strict_power_target),
      power_at_max_N = power_df$power[which.max(power_df$n_subjects)],
      stringsAsFactors = FALSE
    )
  } else {
    high_label <- scenario_df$scenario_label[nrow(scenario_df)]
    sub <- power_df[as.character(power_df$scenario_label) == as.character(high_label), , drop = FALSE]
    summary_row <- data.frame(
      hypothesis = config$hypothesis,
      model = config$label,
      N_min_for_90 = n_for_target_power(sub, strict_power_target),
      power_at_max_N = sub$power[which.max(sub$n_subjects)],
      stringsAsFactors = FALSE
    )
  }

  if (save_figures && nrow(scenario_df) > 1L) {
    message(sprintf("  Saving figures: %s", config$output_prefix))
    save_power_figures(
      power_df = power_df,
      scenario_levels = unique(scenario_df$scenario_display),
      subject_breaks = subject_breaks,
      strict_power_target = strict_power_target,
      figure_output_dir = figure_output_dir,
      output_prefix = config$output_prefix,
      plot_title = config$plot_title,
      y_legend_title = "Residuals SD"
    )
  }

  list(power_df = power_df, summary_row = summary_row)
}
