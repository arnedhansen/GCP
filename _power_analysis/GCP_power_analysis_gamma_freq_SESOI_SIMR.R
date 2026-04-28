## Standalone GCP simulation-based power analysis for gamma frequency (SESOI-only).

ensure_packages <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    user_lib <- Sys.getenv("R_LIBS_USER")
    if (!nzchar(user_lib)) {
      user_lib <- file.path(
        path.expand("~"),
        "R",
        paste0(R.version$platform, "-library"),
        paste(R.version$major, strsplit(R.version$minor, ".", fixed = TRUE)[[1]][1], sep = ".")
      )
    }
    if (!dir.exists(user_lib)) {
      dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
    }
    .libPaths(unique(c(user_lib, .libPaths())))
    install.packages(missing, repos = "https://cloud.r-project.org", lib = user_lib)
  }
}

ensure_packages(c("lme4"))
suppressPackageStartupMessages(library(lme4))

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
    list(gcp_root = file.path("W:/Students/Arne/GCP"))
  } else {
    list(gcp_root = file.path("/Volumes/g_psyplafor_methlab$/Students/Arne/GCP"))
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
  stop("power_parameter_manifest.csv not found. Run _power_analysis/GCP_power_analysis_pilot_stats.R first.")
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
  if (is.finite(med_trials)) as.integer(round(med_trials)) else default_value
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
  data.frame(Subject = Subject, contrast_num_c = contrast_num_c)
}

simulate_outcome <- function(dat, scenario, sim_col, target_term) {
  n_subjects <- nlevels(dat$Subject)
  random_intercepts <- rnorm(n_subjects, mean = 0, sd = scenario$random_intercept_sd * scenario$random_intercept_sd_multiplier)
  random_slopes <- rnorm(n_subjects, mean = 0, sd = scenario$random_slope_sd * scenario$random_slope_sd_multiplier)
  x <- dat$contrast_num_c
  mu <- scenario$outcome_mean + random_intercepts[dat$Subject] + random_slopes[dat$Subject] * x
  if (target_term == "contrast_num_c") {
    mu <- mu + scenario$beta_raw * x
  } else {
    stop("Unsupported target term for standalone SESOI script: ", target_term)
  }
  y <- mu + rnorm(nrow(dat), mean = 0, sd = scenario$residual_sd * scenario$residual_sd_multiplier)
  dat[[sim_col]] <- y
  dat
}

fit_model_and_extract <- function(dat, model_formula, target_term, alpha) {
  fit <- suppressMessages(lmer(model_formula, data = dat, REML = FALSE))
  cf <- as.data.frame(summary(fit)$coefficients)
  cf$term <- rownames(cf)
  target_row <- cf[cf$term == target_term, , drop = FALSE]
  if (nrow(target_row) != 1) {
    stop("Target term not found in model coefficients: ", target_term)
  }
  if ("Pr(>|t|)" %in% names(target_row)) {
    p_value <- as.numeric(target_row[["Pr(>|t|)"]])
  } else if ("t value" %in% names(target_row)) {
    p_value <- 2 * stats::pnorm(abs(as.numeric(target_row[["t value"]])), lower.tail = FALSE)
  } else {
    stop("No p-value or t-value column found in model summary.")
  }
  is.finite(p_value) && p_value < alpha
}

estimate_power <- function(scenario, n_subjects, nsim, trials_per_condition, sim_col, target_term, model_formula, alpha, verbose = TRUE) {
  log_progress("Start SESOI | N=", n_subjects, " | nsim=", nsim, verbose = verbose)
  rejects <- logical(nsim)
  for (i in seq_len(nsim)) {
    dat <- make_subject_condition_design(n_subjects, trials_per_condition)
    dat <- simulate_outcome(dat, scenario, sim_col = sim_col, target_term = target_term)
    rejects[i] <- fit_model_and_extract(dat, model_formula = model_formula, target_term = target_term, alpha = alpha)
  }
  power <- mean(rejects)
  se <- sqrt(power * (1 - power) / nsim)
  out <- data.frame(
    scenario_label = scenario$scenario_label,
    scenario_role = scenario$scenario_role,
    n_subjects = n_subjects,
    power = power,
    lower = pmax(0, power - 1.96 * se),
    upper = pmin(1, power + 1.96 * se),
    nsim = nsim
  )
  log_progress("Done SESOI | N=", n_subjects, " | power=", sprintf("%.3f", out$power), verbose = verbose)
  out
}

run_sesoi_only <- function() {
  cfg <- list(
    manifest_outcome = "PeakFrequency",
    sim_col = "gamma_frequency",
    target_term = "contrast_num_c",
    model_formula = gamma_frequency ~ contrast_num_c + (1 + contrast_num_c | Subject),
    alpha = 0.05,
    nsim = 1000,
    subject_breaks = seq(20, 140, by = 10),
    trials_per_condition = resolve_trials_per_condition(default_value = 160),
    seed = 123,
    strict_power_target = 0.90,
    file_prefix = "GCP_power_analysis_gamma_frequency_SESOI_only",
    verbose = TRUE
  )

  set.seed(cfg$seed)
  output_dir <- resolve_power_output_dir()
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  manifest <- read.csv(resolve_manifest_path(output_dir), stringsAsFactors = FALSE)
  scenario <- manifest[
    manifest$outcome == cfg$manifest_outcome & manifest$scenario_role == "sesoi",
    ,
    drop = FALSE
  ]
  if (nrow(scenario) != 1) {
    stop("Expected exactly one SESOI scenario row for ", cfg$manifest_outcome, ".")
  }

  rows <- lapply(cfg$subject_breaks, function(n_subjects) {
    estimate_power(
      scenario = scenario,
      n_subjects = n_subjects,
      nsim = cfg$nsim,
      trials_per_condition = cfg$trials_per_condition,
      sim_col = cfg$sim_col,
      target_term = cfg$target_term,
      model_formula = cfg$model_formula,
      alpha = cfg$alpha,
      verbose = cfg$verbose
    )
  })
  power_df <- do.call(rbind, rows)
  power_df$meets_target_90 <- power_df$power >= cfg$strict_power_target

  write.csv(power_df, file.path(output_dir, paste0(cfg$file_prefix, "_curve.csv")), row.names = FALSE)
  write.csv(
    data.frame(
      outcome = cfg$manifest_outcome,
      scenario_role = "sesoi",
      strict_power_target = cfg$strict_power_target,
      alpha = cfg$alpha,
      nsim = cfg$nsim,
      stringsAsFactors = FALSE
    ),
    file.path(output_dir, paste0(cfg$file_prefix, "_assumptions_manifest.csv")),
    row.names = FALSE
  )
  power_df
}

result <- run_sesoi_only()
print(result)
