## Standalone parallel-only GCP simulation-based power analysis for gamma frequency (SESOI-only).

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
    list(gcp_root = file.path("W:/Students/Arne/GCP"))
  } else {
    list(gcp_root = file.path("/Volumes/g_psyplafor_methlab$/Students/Arne/GCP"))
  }
}

resolve_power_output_dir <- function() {
  roots <- resolve_project_roots()
  preferred <- file.path(roots$gcp_root, "figures", "power_analysis")
  fallback <- file.path(getwd(), "_power_analysis", "outputs")
  if (dir.exists(dirname(preferred))) preferred else fallback
}

resolve_trials_per_condition <- function(default_value = 160) {
  roots <- resolve_project_roots()
  candidates <- c(
    file.path(roots$gcp_root, "data", "features", "GCP_eeg_GED_gamma_metrics_trials.csv"),
    file.path(getwd(), "data", "features", "GCP_eeg_GED_gamma_metrics_trials.csv")
  )
  input_file <- candidates[file.exists(candidates)][1]
  dat <- read.csv(input_file, stringsAsFactors = FALSE)
  counts <- aggregate(rep(1, nrow(dat)) ~ Subject + Contrast, data = dat, FUN = length)
  names(counts)[3] <- "n_trials"
  as.integer(round(stats::median(counts$n_trials, na.rm = TRUE)))
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

simulate_outcome <- function(dat, scenario, sim_col) {
  as_num <- function(x) suppressWarnings(as.numeric(x[1]))
  n_subjects <- nlevels(dat$Subject)
  ri_sd <- as_num(scenario$random_intercept_sd)
  rs_sd <- as_num(scenario$random_slope_sd)
  e_sd <- as_num(scenario$residual_sd)
  beta_raw <- as_num(scenario$beta_raw)
  mu0 <- as_num(scenario$outcome_mean)
  ri_mult <- as_num(scenario$random_intercept_sd_multiplier)
  rs_mult <- as_num(scenario$random_slope_sd_multiplier)
  e_mult <- as_num(scenario$residual_sd_multiplier)
  if (!is.finite(ri_mult)) ri_mult <- 1
  if (!is.finite(rs_mult)) rs_mult <- 1
  if (!is.finite(e_mult)) e_mult <- 1
  random_intercepts <- rnorm(n_subjects, mean = 0, sd = ri_sd * ri_mult)
  random_slopes <- rnorm(n_subjects, mean = 0, sd = rs_sd * rs_mult)
  x <- dat$contrast_num_c
  mu <- mu0 + random_intercepts[dat$Subject] + random_slopes[dat$Subject] * x + beta_raw * x
  y <- mu + rnorm(nrow(dat), mean = 0, sd = e_sd * e_mult)
  dat[[sim_col]] <- y
  dat
}

fit_model_and_extract <- function(dat, model_formula, target_term, alpha) {
  fit <- suppressMessages(lmer(model_formula, data = dat, REML = FALSE))
  cf <- as.data.frame(summary(fit)$coefficients)
  cf$term <- rownames(cf)
  target_row <- cf[cf$term == target_term, , drop = FALSE]
  if ("Pr(>|t|)" %in% names(target_row)) {
    p_value <- as.numeric(target_row[["Pr(>|t|)"]])
  } else {
    p_value <- 2 * stats::pnorm(abs(as.numeric(target_row[["t value"]])), lower.tail = FALSE)
  }
  is.finite(p_value) && p_value < alpha
}

estimate_power_chunk <- function(
    chunk_nsim,
    scenario,
    n_subjects,
    trials_per_condition,
    sim_col,
    model_formula,
    target_term,
    alpha) {
  rejects <- 0L
  for (i in seq_len(chunk_nsim)) {
    dat <- make_subject_condition_design(n_subjects, trials_per_condition)
    dat <- simulate_outcome(dat, scenario, sim_col = sim_col)
    rejects <- rejects + as.integer(fit_model_and_extract(dat, model_formula, target_term, alpha))
  }
  list(chunk_nsim = as.integer(chunk_nsim), rejects = as.integer(rejects))
}

estimate_power_parallel <- function(
    scenario,
    n_subjects,
    nsim,
    trials_per_condition,
    sim_col,
    model_formula,
    target_term,
    alpha,
    cl,
    verbose = TRUE,
    round_chunk_nsim = 1L) {
  n_workers <- length(cl)
  log_progress("Start SESOI parallel | N=", n_subjects, " | nsim=", nsim, " | workers=", n_workers, verbose = verbose)
  total_nsim <- 0L
  total_rejects <- 0L
  remaining <- as.integer(nsim)
  round_idx <- 0L

  while (remaining > 0L) {
    round_idx <- round_idx + 1L
    workers_this_round <- min(n_workers, remaining)
    chunk_sizes <- rep(round_chunk_nsim, workers_this_round)
    max_budget <- workers_this_round * round_chunk_nsim
    if (max_budget > remaining) {
      overflow <- max_budget - remaining
      for (k in seq_len(overflow)) {
        idx <- ((k - 1L) %% workers_this_round) + 1L
        chunk_sizes[idx] <- chunk_sizes[idx] - 1L
      }
    }
    chunk_sizes <- chunk_sizes[chunk_sizes > 0L]

    chunk_results <- parallel::parLapply(
      cl,
      chunk_sizes,
      fun = function(chunk_nsim, scenario, n_subjects, trials_per_condition, sim_col, model_formula, target_term, alpha) {
        estimate_power_chunk(
          chunk_nsim = chunk_nsim,
          scenario = scenario,
          n_subjects = n_subjects,
          trials_per_condition = trials_per_condition,
          sim_col = sim_col,
          model_formula = model_formula,
          target_term = target_term,
          alpha = alpha
        )
      },
      scenario = scenario,
      n_subjects = n_subjects,
      trials_per_condition = trials_per_condition,
      sim_col = sim_col,
      model_formula = model_formula,
      target_term = target_term,
      alpha = alpha
    )

    round_n <- sum(vapply(chunk_results, function(x) as.integer(x$chunk_nsim), integer(1)))
    total_nsim <- total_nsim + round_n
    total_rejects <- total_rejects + sum(vapply(chunk_results, function(x) as.integer(x$rejects), integer(1)))
    remaining <- nsim - total_nsim
    interim_power <- total_rejects / total_nsim
    log_progress(
      "Heartbeat SESOI parallel | N=", n_subjects,
      " | round=", round_idx,
      " | processed=", total_nsim, "/", nsim,
      " | interim_power=", sprintf("%.3f", interim_power),
      verbose = verbose
    )
  }

  power <- total_rejects / total_nsim
  se <- sqrt(power * (1 - power) / total_nsim)
  out <- data.frame(
    scenario_label = scenario$scenario_label,
    scenario_role = scenario$scenario_role,
    n_subjects = n_subjects,
    power = power,
    lower = pmax(0, power - 1.96 * se),
    upper = pmin(1, power + 1.96 * se),
    nsim = total_nsim
  )
  log_progress("Done SESOI parallel | N=", n_subjects, " | power=", sprintf("%.3f", out$power), verbose = verbose)
  out
}

run_sesoi_only <- function() {
  cfg <- list(
    sim_col = "gamma_frequency",
    target_term = "contrast_num_c",
    model_formula = gamma_frequency ~ contrast_num_c + (1 + contrast_num_c | Subject),
    alpha = 0.05,
    nsim = 1000,
    subject_breaks = seq(20, 60, by = 5),
    trials_per_condition = resolve_trials_per_condition(default_value = 160),
    seed = 123,
    strict_power_target = 0.90,
    parallel_workers = 8,
    parallel_round_chunk_nsim = 1,
    file_prefix = "GCP_power_analysis_gamma_frequency_SESOI_only",
    plot_title = "Power Curve: Gamma Frequency SESOI (contrast_num_c)",
    verbose = TRUE
  )

  set.seed(cfg$seed)
  output_dir <- resolve_power_output_dir()
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  scenario <- data.frame(
    scenario_label = "SESOI",
    scenario_role = "sesoi",
    beta_raw = 0.05,
    outcome_mean = 0.00,
    random_intercept_sd = 0.20,
    random_slope_sd = 0.10,
    residual_sd = 1.00,
    random_intercept_sd_multiplier = 1.00,
    random_slope_sd_multiplier = 1.00,
    residual_sd_multiplier = 1.00,
    stringsAsFactors = FALSE
  )

  cl <- parallel::makeCluster(as.integer(cfg$parallel_workers), type = "PSOCK")
  on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::clusterEvalQ(cl, suppressPackageStartupMessages(library(lme4)))
  parallel::clusterSetRNGStream(cl, iseed = cfg$seed)
  parallel::clusterExport(
    cl,
    varlist = c(
      "make_subject_condition_design",
      "simulate_outcome",
      "fit_model_and_extract",
      "estimate_power_chunk"
    ),
    envir = environment()
  )

  rows <- lapply(cfg$subject_breaks, function(n_subjects) {
    estimate_power_parallel(
      scenario = scenario,
      n_subjects = n_subjects,
      nsim = cfg$nsim,
      trials_per_condition = cfg$trials_per_condition,
      sim_col = cfg$sim_col,
      model_formula = cfg$model_formula,
      target_term = cfg$target_term,
      alpha = cfg$alpha,
      cl = cl,
      verbose = cfg$verbose,
      round_chunk_nsim = cfg$parallel_round_chunk_nsim
    )
  })

  power_df <- do.call(rbind, rows)
  power_df$meets_target_90 <- power_df$power >= cfg$strict_power_target
  write.csv(power_df, file.path(output_dir, paste0(cfg$file_prefix, "_curve.csv")), row.names = FALSE)

  x_breaks <- sort(unique(cfg$subject_breaks))
  plot_obj <- ggplot(power_df, aes(x = .data$n_subjects, y = .data$power)) +
    geom_point(size = 2, color = "blue") +
    geom_line(color = "blue", linetype = "dotted") +
    geom_errorbar(aes(ymin = .data$lower, ymax = .data$upper), width = 3, color = "blue", alpha = 0.5) +
    geom_hline(yintercept = cfg$strict_power_target, linetype = "dashed", color = "grey", linewidth = 0.5) +
    scale_x_continuous(breaks = x_breaks) +
    scale_y_continuous(
      labels = scales::percent_format(),
      limits = c(0, 1),
      breaks = c(0, 0.25, 0.50, 0.75, 0.90, 1)
    ) +
    labs(
      x = "Subjects",
      y = "Power",
      title = cfg$plot_title
    ) +
    theme_minimal(base_size = 15) +
    theme(
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
      axis.line = element_line(colour = "black", linewidth = 1),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.ticks.length = grid::unit(0.2, "cm"),
      axis.ticks = element_line(linewidth = 0.5),
      axis.ticks.x = element_line(colour = "black", linewidth = 0.5, lineend = "square"),
      axis.ticks.y = element_line(colour = "black", linewidth = 0.5, lineend = "square"),
      axis.text.x = element_text(margin = margin(t = 10)),
      axis.text.y = element_text(margin = margin(r = 10))
    )

  png(
    file = file.path(output_dir, paste0(cfg$file_prefix, ".png")),
    width = 2200, height = 1400, res = 300
  )
  print(plot_obj)
  dev.off()

  power_df
}

result <- run_sesoi_only()
print(result)
