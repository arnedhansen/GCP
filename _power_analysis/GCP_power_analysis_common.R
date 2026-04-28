## Shared utilities for GCP power analysis scripts.

ensure_packages <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    user_lib <- Sys.getenv("R_LIBS_USER")
    if (!nzchar(user_lib)) {
      user_lib <- file.path(path.expand("~"), "R", paste0(R.version$platform, "-library"), paste(R.version$major, strsplit(R.version$minor, ".", fixed = TRUE)[[1]][1], sep = "."))
    }
    if (!dir.exists(user_lib)) {
      dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
    }
    .libPaths(unique(c(user_lib, .libPaths())))
    install.packages(missing, repos = "https://cloud.r-project.org", lib = user_lib)
    still_missing <- missing[!vapply(missing, requireNamespace, logical(1), quietly = TRUE)]
    if (length(still_missing) > 0) {
      stop(
        "Failed to install required packages: ",
        paste(still_missing, collapse = ", "),
        ". Try manually in R: install.packages(c(",
        paste(sprintf("'%s'", still_missing), collapse = ", "),
        "), repos='https://cloud.r-project.org')"
      )
    }
  }
}

ensure_packages(c("lme4", "ggplot2", "scales"))
suppressPackageStartupMessages({
  library(lme4)
  library(ggplot2)
})

log_progress <- function(..., verbose = TRUE) {
  if (!isTRUE(verbose)) {
    return(invisible(NULL))
  }
  cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = "")
  flush.console()
  invisible(NULL)
}

resolve_project_roots <- function() {
  if (.Platform$OS.type == "windows") {
    list(
      gcp_root = file.path("W:/Students/Arne/GCP"),
      share_root = "W:/"
    )
  } else {
    list(
      gcp_root = file.path("/Volumes/g_psyplafor_methlab$/Students/Arne/GCP"),
      share_root = "/Volumes/g_psyplafor_methlab$/"
    )
  }
}

resolve_power_output_dir <- function() {
  roots <- resolve_project_roots()
  preferred <- file.path(roots$gcp_root, "figures", "power_analysis")
  fallback <- file.path(getwd(), "_power_analysis", "outputs")
  if (dir.exists(dirname(preferred))) {
    return(preferred)
  }
  fallback
}

resolve_manifest_path <- function(output_dir = resolve_power_output_dir()) {
  candidate_a <- file.path(output_dir, "pilot_stats", "power_parameter_manifest.csv")
  candidate_b <- file.path(getwd(), "_power_analysis", "outputs", "pilot_stats", "power_parameter_manifest.csv")
  if (file.exists(candidate_a)) {
    return(candidate_a)
  }
  if (file.exists(candidate_b)) {
    return(candidate_b)
  }
  stop(
    "power_parameter_manifest.csv not found. Run ",
    "_power_analysis/GCP_power_analysis_pilot_stats.R first."
  )
}

resolve_trials_per_condition <- function(default_value = 160) {
  roots <- resolve_project_roots()
  candidates <- c(
    file.path(roots$gcp_root, "data", "features", "GCP_eeg_GED_gamma_metrics_trials.csv"),
    file.path(getwd(), "data", "features", "GCP_eeg_GED_gamma_metrics_trials.csv")
  )
  input_file <- candidates[file.exists(candidates)][1]
  if (!file.exists(input_file)) {
    return(default_value)
  }
  dat <- read.csv(input_file, stringsAsFactors = FALSE)
  if (!all(c("Subject", "Contrast") %in% names(dat))) {
    return(default_value)
  }
  counts <- aggregate(rep(1, nrow(dat)) ~ Subject + Contrast, data = dat, FUN = length)
  names(counts)[3] <- "n_trials"
  med_trials <- stats::median(counts$n_trials, na.rm = TRUE)
  if (is.finite(med_trials)) {
    return(as.integer(round(med_trials)))
  }
  default_value
}

make_subject_condition_design <- function(n_subjects, trials_per_condition, contrast_levels = c("25", "50", "75", "100")) {
  Subject <- factor(rep(seq_len(n_subjects), each = length(contrast_levels) * trials_per_condition))
  contrast <- factor(
    rep(rep(contrast_levels, each = trials_per_condition), times = n_subjects),
    levels = contrast_levels,
    ordered = TRUE
  )
  contrast_num <- as.numeric(as.character(contrast))
  contrast_num_c <- as.numeric(scale(contrast_num, center = TRUE, scale = TRUE))
  contrast_num_c2 <- contrast_num_c^2
  data.frame(
    Subject = Subject,
    contrast = contrast,
    contrast_num = contrast_num,
    contrast_num_c = contrast_num_c,
    contrast_num_c2 = contrast_num_c2
  )
}

simulate_outcome <- function(dat, scenario, outcome_cfg) {
  n_subjects <- nlevels(dat$Subject)
  random_intercepts <- rnorm(n_subjects, mean = 0, sd = scenario$random_intercept_sd)
  random_slopes <- rnorm(n_subjects, mean = 0, sd = scenario$random_slope_sd)
  x <- dat$contrast_num_c
  mu <- scenario$outcome_mean +
    random_intercepts[dat$Subject] +
    random_slopes[dat$Subject] * x
  if (outcome_cfg$target_term == "contrast_num_c") {
    mu <- mu + scenario$beta_raw * x
  } else if (outcome_cfg$target_term == "contrast_num_c2") {
    mu <- mu + outcome_cfg$linear_nuisance_beta * x + scenario$beta_raw * (x^2)
  } else {
    stop("Unsupported target term: ", outcome_cfg$target_term)
  }
  y <- mu + rnorm(nrow(dat), mean = 0, sd = scenario$residual_sd)
  dat[[outcome_cfg$sim_col]] <- y
  dat
}

fit_model_and_extract <- function(dat, outcome_cfg, alpha) {
  form <- outcome_cfg$model_formula
  warn_bucket <- character(0)
  fit <- withCallingHandlers(
    suppressMessages(lmer(form, data = dat, REML = FALSE)),
    warning = function(w) {
      warn_bucket <<- c(warn_bucket, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  cf <- as.data.frame(summary(fit)$coefficients)
  cf$term <- rownames(cf)
  p_col <- grep("^Pr\\(", names(cf), value = TRUE)
  target_row <- cf[cf$term == outcome_cfg$target_term, , drop = FALSE]
  if (nrow(target_row) != 1) {
    stop("Target term not found in model coefficients: ", outcome_cfg$target_term)
  }
  if (length(p_col) > 0) {
    p_value <- as.numeric(target_row[[p_col[1]]])
  } else if ("t value" %in% names(target_row)) {
    ## lme4 does not report p-values by default; use normal approximation.
    p_value <- 2 * stats::pnorm(abs(as.numeric(target_row[["t value"]])), lower.tail = FALSE)
  } else {
    stop("No p-value or t-value column found in model summary.")
  }
  list(
    p_value = p_value,
    reject = is.finite(p_value) && p_value < alpha,
    singular = isSingular(fit, tol = 1e-4),
    warnings = warn_bucket
  )
}

estimate_power <- function(
    scenario,
    outcome_cfg,
    n_subjects,
    nsim,
    trials_per_condition,
    alpha,
    verbose = TRUE,
    progress_every = NULL,
    context_label = NULL) {
  rejects <- logical(nsim)
  singular_hits <- logical(nsim)
  warning_count <- integer(nsim)

  if (is.null(progress_every)) {
    progress_every <- max(1L, as.integer(floor(nsim / 10L)))
  }
  if (is.null(context_label) || !nzchar(context_label)) {
    context_label <- sprintf("scenario=%s n=%d", scenario$scenario_label, n_subjects)
  }
  log_progress("Start: ", context_label, " | nsim=", nsim, verbose = verbose)

  for (i in seq_len(nsim)) {
    dat <- make_subject_condition_design(n_subjects, trials_per_condition)
    dat <- simulate_outcome(dat, scenario, outcome_cfg)
    fit_out <- fit_model_and_extract(dat, outcome_cfg, alpha)
    rejects[i] <- fit_out$reject
    singular_hits[i] <- fit_out$singular
    warning_count[i] <- length(fit_out$warnings)

    if (i %% progress_every == 0L || i == nsim) {
      interim_power <- mean(rejects[seq_len(i)])
      log_progress(
        "Progress: ", context_label, " | ",
        i, "/", nsim,
        " sims | interim_power=", sprintf("%.3f", interim_power),
        verbose = verbose
      )
    }
  }

  power <- mean(rejects)
  se <- sqrt(power * (1 - power) / nsim)
  out <- data.frame(
    scenario_label = scenario$scenario_label,
    scenario_role = scenario$scenario_role,
    effect_source = scenario$effect_source,
    target_term = outcome_cfg$target_term,
    n_subjects = n_subjects,
    power = power,
    lower = pmax(0, power - 1.96 * se),
    upper = pmin(1, power + 1.96 * se),
    mc_se = se,
    nsim = nsim,
    singular_rate = mean(singular_hits),
    warning_rate = mean(warning_count > 0),
    warnings_per_fit = mean(warning_count)
  )
  log_progress(
    "Done: ", context_label,
    " | power=", sprintf("%.3f", out$power),
    " [", sprintf("%.3f", out$lower), ", ", sprintf("%.3f", out$upper), "]",
    verbose = verbose
  )
  out
}

run_validation_suite <- function(base_scenario, outcome_cfg, cfg) {
  verbose <- isTRUE(cfg$verbose)
  log_progress("Validation suite started.", verbose = verbose)
  null_scenario <- base_scenario
  null_scenario$scenario_label <- "Validation_Null"
  null_scenario$scenario_role <- "validation"
  null_scenario$effect_source <- "forced_null"
  null_scenario$beta_raw <- 0

  null_res <- estimate_power(
    scenario = null_scenario,
    outcome_cfg = outcome_cfg,
    n_subjects = cfg$validation_n,
    nsim = cfg$validation_nsim,
    trials_per_condition = cfg$trials_per_condition,
    alpha = cfg$alpha,
    verbose = verbose,
    progress_every = cfg$progress_every,
    context_label = "validation:null_type_I"
  )

  large_n_res <- estimate_power(
    scenario = base_scenario,
    outcome_cfg = outcome_cfg,
    n_subjects = max(cfg$subject_breaks),
    nsim = cfg$validation_nsim,
    trials_per_condition = cfg$trials_per_condition,
    alpha = cfg$alpha,
    verbose = verbose,
    progress_every = cfg$progress_every,
    context_label = "validation:large_N"
  )
  large_n_res$scenario_label <- "Validation_LargeN"
  large_n_res$scenario_role <- "validation"

  monotonic_n <- cfg$subject_breaks
  monotonic_power <- numeric(length(monotonic_n))
  for (i in seq_along(monotonic_n)) {
    r <- estimate_power(
      scenario = base_scenario,
      outcome_cfg = outcome_cfg,
      n_subjects = monotonic_n[i],
      nsim = cfg$validation_nsim,
      trials_per_condition = cfg$trials_per_condition,
      alpha = cfg$alpha,
      verbose = FALSE,
      progress_every = cfg$progress_every
    )
    monotonic_power[i] <- r$power
  }
  monotonic_pass <- all(diff(monotonic_power) > -0.05)
  monotonicity <- data.frame(
    check = "monotonicity",
    passed = monotonic_pass,
    details = paste(round(monotonic_power, 3), collapse = " -> "),
    stringsAsFactors = FALSE
  )

  checks <- rbind(
    data.frame(
      check = "null_type_I_rate",
      passed = abs(null_res$power - cfg$alpha) <= cfg$alpha_tolerance,
      details = sprintf("observed=%.3f expected~%.2f", null_res$power, cfg$alpha),
      stringsAsFactors = FALSE
    ),
    data.frame(
      check = "large_n_power",
      passed = large_n_res$power >= cfg$large_n_min_power,
      details = sprintf("observed=%.3f threshold=%.2f", large_n_res$power, cfg$large_n_min_power),
      stringsAsFactors = FALSE
    ),
    monotonicity
  )

  list(
    validation_curve = rbind(null_res, large_n_res),
    checks = checks
  )
}

compute_recommendations <- function(power_df, strict_base = 0.90, strict_pessimistic = 0.80, feasible_cap = 80) {
  base <- power_df[power_df$scenario_role == "base", c("n_subjects", "power")]
  pessimistic <- power_df[power_df$scenario_role == "pessimistic", c("n_subjects", "power")]
  names(base)[2] <- "power_base"
  names(pessimistic)[2] <- "power_pess"
  merged <- merge(base, pessimistic, by = "n_subjects", all = FALSE)
  merged <- merged[order(merged$n_subjects), ]

  strict_hits <- merged[merged$power_base >= strict_base & merged$power_pess >= strict_pessimistic, , drop = FALSE]
  if (nrow(strict_hits) == 0) {
    strict_row <- data.frame(
      recommendation_type = "strict",
      n_subjects = NA_real_,
      power_base = NA_real_,
      power_pess = NA_real_,
      shortfall_base = NA_real_,
      shortfall_pess = NA_real_,
      note = "No tested sample size met strict criteria."
    )
  } else {
    strict_row <- strict_hits[1, , drop = FALSE]
    strict_row$recommendation_type <- "strict"
    strict_row$shortfall_base <- pmax(0, strict_base - strict_row$power_base)
    strict_row$shortfall_pess <- pmax(0, strict_pessimistic - strict_row$power_pess)
    strict_row$note <- "Meets base>=0.90 and pessimistic>=0.80."
    strict_row <- strict_row[, c(
      "recommendation_type", "n_subjects", "power_base", "power_pess",
      "shortfall_base", "shortfall_pess", "note"
    )]
  }

  feasible_pool <- merged[merged$n_subjects <= feasible_cap, , drop = FALSE]
  if (nrow(feasible_pool) == 0) {
    feasible_row <- data.frame(
      recommendation_type = "feasible",
      n_subjects = NA_real_,
      power_base = NA_real_,
      power_pess = NA_real_,
      shortfall_base = NA_real_,
      shortfall_pess = NA_real_,
      note = "No sample size falls under feasibility cap."
    )
  } else {
    feasible_pool$criterion_distance <- pmin(
      feasible_pool$power_base / strict_base,
      feasible_pool$power_pess / strict_pessimistic
    )
    feasible_pool <- feasible_pool[order(-feasible_pool$criterion_distance, feasible_pool$n_subjects), ]
    feasible_row <- feasible_pool[1, , drop = FALSE]
    feasible_row$recommendation_type <- "feasible"
    feasible_row$shortfall_base <- pmax(0, strict_base - feasible_row$power_base)
    feasible_row$shortfall_pess <- pmax(0, strict_pessimistic - feasible_row$power_pess)
    feasible_row$note <- sprintf("Best design under feasibility cap n<=%d.", feasible_cap)
    feasible_row <- feasible_row[, c(
      "recommendation_type", "n_subjects", "power_base", "power_pess",
      "shortfall_base", "shortfall_pess", "note"
    )]
  }

  rbind(strict_row, feasible_row)
}

run_outcome_power_analysis <- function(outcome_cfg) {
  verbose <- if (!is.null(outcome_cfg$verbose)) isTRUE(outcome_cfg$verbose) else TRUE
  progress_every <- if (!is.null(outcome_cfg$progress_every)) as.integer(outcome_cfg$progress_every) else NULL

  set.seed(outcome_cfg$seed)

  run_mode <- outcome_cfg$run_mode
  nsim <- if (identical(run_mode, "reporting")) outcome_cfg$nsim_reporting else outcome_cfg$nsim_exploratory
  log_progress(
    "Run started: outcome=", outcome_cfg$manifest_outcome,
    " | mode=", run_mode,
    " | nsim=", nsim,
    " | seed=", outcome_cfg$seed,
    verbose = verbose
  )

  output_dir <- resolve_power_output_dir()
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  manifest <- read.csv(resolve_manifest_path(output_dir), stringsAsFactors = FALSE)

  scenario_manifest <- manifest[manifest$outcome == outcome_cfg$manifest_outcome, , drop = FALSE]
  if (nrow(scenario_manifest) == 0) {
    stop("No manifest rows found for outcome ", outcome_cfg$manifest_outcome)
  }

  subject_breaks <- outcome_cfg$subject_breaks
  all_rows <- list()
  row_idx <- 1L
  for (i in seq_len(nrow(scenario_manifest))) {
    sc <- scenario_manifest[i, , drop = FALSE]
    log_progress(
      "Scenario ", i, "/", nrow(scenario_manifest),
      ": label=", sc$scenario_label,
      " role=", sc$scenario_role,
      " source=", sc$effect_source,
      verbose = verbose
    )
    for (n in subject_breaks) {
      all_rows[[row_idx]] <- estimate_power(
        scenario = sc,
        outcome_cfg = outcome_cfg,
        n_subjects = n,
        nsim = nsim,
        trials_per_condition = outcome_cfg$trials_per_condition,
        alpha = outcome_cfg$alpha,
        verbose = verbose,
        progress_every = progress_every,
        context_label = sprintf("%s | N=%d", sc$scenario_label, n)
      )
      row_idx <- row_idx + 1L
    }
  }
  power_df <- do.call(rbind, all_rows)

  base_sc <- scenario_manifest[scenario_manifest$scenario_role == "base", , drop = FALSE]
  if (nrow(base_sc) != 1) {
    stop("Expected exactly one base scenario for ", outcome_cfg$manifest_outcome)
  }
  validation <- run_validation_suite(base_sc, outcome_cfg, list(
    validation_n = outcome_cfg$validation_n,
    validation_nsim = outcome_cfg$validation_nsim,
    trials_per_condition = outcome_cfg$trials_per_condition,
    alpha = outcome_cfg$alpha,
    alpha_tolerance = outcome_cfg$alpha_tolerance,
    large_n_min_power = outcome_cfg$large_n_min_power,
    subject_breaks = subject_breaks,
    verbose = verbose,
    progress_every = progress_every
  ))

  recommendations <- compute_recommendations(
    power_df = power_df,
    strict_base = outcome_cfg$strict_power_base,
    strict_pessimistic = outcome_cfg$strict_power_pessimistic,
    feasible_cap = outcome_cfg$feasible_subject_cap
  )

  assumptions_manifest <- data.frame(
    outcome = outcome_cfg$manifest_outcome,
    target_term = outcome_cfg$target_term,
    alpha = outcome_cfg$alpha,
    run_mode = run_mode,
    nsim = nsim,
    nsim_validation = outcome_cfg$validation_nsim,
    trials_per_condition = outcome_cfg$trials_per_condition,
    strict_power_base = outcome_cfg$strict_power_base,
    strict_power_pessimistic = outcome_cfg$strict_power_pessimistic,
    feasible_subject_cap = outcome_cfg$feasible_subject_cap,
    seed = outcome_cfg$seed,
    subject_breaks = paste(subject_breaks, collapse = ","),
    stringsAsFactors = FALSE
  )

  file_prefix <- outcome_cfg$file_prefix
  write.csv(power_df, file.path(output_dir, paste0(file_prefix, "_curve.csv")), row.names = FALSE)
  write.csv(validation$validation_curve, file.path(output_dir, paste0(file_prefix, "_validation_curve.csv")), row.names = FALSE)
  write.csv(validation$checks, file.path(output_dir, paste0(file_prefix, "_validation_checks.csv")), row.names = FALSE)
  write.csv(recommendations, file.path(output_dir, paste0(file_prefix, "_recommendations.csv")), row.names = FALSE)
  write.csv(assumptions_manifest, file.path(output_dir, paste0(file_prefix, "_assumptions_manifest.csv")), row.names = FALSE)

  plot_obj <- ggplot(power_df, aes(
    x = n_subjects, y = power,
    color = scenario_label, linetype = scenario_role
  )) +
    geom_line(linewidth = 1) +
    geom_point(size = 1.8) +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 1.5, alpha = 0.2) +
    geom_hline(yintercept = outcome_cfg$strict_power_base, linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = outcome_cfg$strict_power_pessimistic, linetype = "dotted", color = "grey50") +
    scale_y_continuous(limits = c(0, 1), labels = scales::percent_format()) +
    labs(
      x = "Subjects",
      y = "Power",
      color = "Scenario",
      linetype = "Scenario role",
      title = outcome_cfg$plot_title
    ) +
    theme_minimal(base_size = 14)

  png(
    file = file.path(output_dir, paste0(file_prefix, ".png")),
    width = 2200, height = 1400, res = 300
  )
  print(plot_obj)
  dev.off()

  saveRDS(
    list(
      scenario_manifest = scenario_manifest,
      power_curve = power_df,
      validation = validation,
      recommendations = recommendations,
      assumptions_manifest = assumptions_manifest
    ),
    file = file.path(output_dir, paste0(file_prefix, "_results.rds"))
  )

  log_progress("Run finished. Outputs written to: ", output_dir, verbose = verbose)
  list(
    power_curve = power_df,
    validation = validation,
    recommendations = recommendations,
    output_dir = output_dir
  )
}
