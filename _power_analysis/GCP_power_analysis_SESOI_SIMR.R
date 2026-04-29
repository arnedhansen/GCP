## Simulates SESOI-based power with mixed models and saves a CSV power curve plus PNG figure

# install.packages("lme4", repos = "https://cloud.r-project.org")
# install.packages("ggplot2", repos = "https://cloud.r-project.org")
# install.packages("scales", repos = "https://cloud.r-project.org")

suppressPackageStartupMessages(library(lme4))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(scales))

resolve_gcp_root <- function() {
  if (.Platform$OS.type == "windows") {
    return("W:/Students/Arne/GCP")
  }
  "/Users/Arne/Documents/GitHub/GCP"
}

runSESOI <- function() {
  # Simulation settings
  seed <- 123
  alpha <- 0.05
  nsim <- 5000
  strict_power_target <- 0.90
  subject_breaks <- seq(5, 75, by = 5)
  contrast_levels <- c("25", "50", "75", "100")
  trials_per_condition <- 160
  sesoi_beta <- 0.10
  outcome_mean <- 0.00
  baseline_random_intercept_sd <- 0.20
  baseline_random_slope_sd <- 0.18
  baseline_residual_sd <- 1.00
  sd_multipliers <- c(0.75, 1.00, 1.25, 1.50)
  trial_missingness_range <- c(0.10, 0.20)
  subject_dropout_range <- c(0.05, 0.10)
  parallel_workers <- 8
  parallel_round_chunk_nsim <- 1
  sim_col <- "gamma_power"
  target_term <- "contrast_num_c"
  model_formula <- gamma_power ~ contrast_num_c + (1 + contrast_num_c | Subject)
  output_prefix <- "GCP_power_analysis_SESOI_linear"
  plot_title <- "Power Analysis: SESOI Linear Slope"
  heat_low_color <- "#BF0D3E"
  heat_high_color <- "#2E7D32"

  gcp_root <- resolve_gcp_root()
  output_dir <- file.path(gcp_root, "figures", "power_analysis")
  set.seed(seed)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # Build sensitivity scenarios by changing one variance block at a time
  make_scenarios <- function() {
    baseline <- data.frame(
      scenario_label = "baseline",
      varied_component = "none",
      multiplier = 1.00,
      random_intercept_sd_value = baseline_random_intercept_sd,
      random_slope_sd_value = baseline_random_slope_sd,
      residual_sd_value = baseline_residual_sd,
      stringsAsFactors = FALSE
    )
    scenario_rows <- list(baseline)
    for (m in sd_multipliers[sd_multipliers != 1]) {
      scenario_rows[[length(scenario_rows) + 1]] <- data.frame(
        scenario_label = sprintf("random_intercept_sd_%.2fx", m),
        varied_component = "random_intercept_sd",
        multiplier = m,
        random_intercept_sd_value = baseline_random_intercept_sd * m,
        random_slope_sd_value = baseline_random_slope_sd,
        residual_sd_value = baseline_residual_sd,
        stringsAsFactors = FALSE
      )
      scenario_rows[[length(scenario_rows) + 1]] <- data.frame(
        scenario_label = sprintf("random_slope_sd_%.2fx", m),
        varied_component = "random_slope_sd",
        multiplier = m,
        random_intercept_sd_value = baseline_random_intercept_sd,
        random_slope_sd_value = baseline_random_slope_sd * m,
        residual_sd_value = baseline_residual_sd,
        stringsAsFactors = FALSE
      )
      scenario_rows[[length(scenario_rows) + 1]] <- data.frame(
        scenario_label = sprintf("residual_sd_%.2fx", m),
        varied_component = "residual_sd",
        multiplier = m,
        random_intercept_sd_value = baseline_random_intercept_sd,
        random_slope_sd_value = baseline_random_slope_sd,
        residual_sd_value = baseline_residual_sd * m,
        stringsAsFactors = FALSE
      )
    }
    do.call(rbind, scenario_rows)
  }

  estimate_power_chunk <- function(
      chunk_nsim,
      n_subjects,
      trials_per_condition,
      sim_col,
      model_formula,
      target_term,
      alpha,
      beta_raw,
      outcome_mean,
      random_intercept_sd,
      random_slope_sd,
      residual_sd,
      trial_missingness_range,
      subject_dropout_range) {
    # Count significant tests for the target fixed effect
    rejects <- 0
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
      mu <- outcome_mean + random_intercepts[dat$Subject] + random_slopes[dat$Subject] * x + beta_raw * x
      dat[[sim_col]] <- mu + rnorm(nrow(dat), mean = 0, sd = residual_sd)

      # Apply trial-level missingness
      trial_missingness_rate <- runif(1, min = trial_missingness_range[1], max = trial_missingness_range[2])
      keep_trial <- stats::runif(nrow(dat)) > trial_missingness_rate
      dat <- dat[keep_trial, , drop = FALSE]
      if (nrow(dat) == 0) {
        next
      }

      # Apply subject-level dropout after trial loss
      subject_dropout_rate <- runif(1, min = subject_dropout_range[1], max = subject_dropout_range[2])
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

      # Skip failed fits without adding a rejection
      fit <- tryCatch(
        suppressMessages(lmer(model_formula, data = dat, REML = FALSE)),
        error = function(e) NULL
      )
      if (is.null(fit)) {
        next
      }
      cf <- as.data.frame(summary(fit)$coefficients)
      cf$term <- rownames(cf)
      target_row <- cf[cf$term == target_term, , drop = FALSE]
      # Fall back to normal-approximation p-value if needed
      if ("Pr(>|t|)" %in% names(target_row)) {
        p_value <- as.numeric(target_row[["Pr(>|t|)"]])
      } else {
        p_value <- 2 * stats::pnorm(abs(as.numeric(target_row[["t value"]])), lower.tail = FALSE)
      }
      rejects <- rejects + as.integer(is.finite(p_value) && p_value < alpha)
    }
    list(chunk_nsim = as.integer(chunk_nsim), rejects = as.integer(rejects))
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
    cat(sprintf("[%s] SCENARIO %s\n", format(Sys.time(), "%H:%M:%S"), scenario$scenario_label))
    flush.console()

    n_rows <- lapply(subject_breaks, function(n_subjects) {
      cat(sprintf("[%s] START scenario=%s N=%d\n", format(Sys.time(), "%H:%M:%S"), scenario$scenario_label, n_subjects))
      flush.console()

      total_nsim <- 0
      total_rejects <- 0
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
          model_formula = model_formula,
          target_term = target_term,
          alpha = alpha,
          beta_raw = sesoi_beta,
          outcome_mean = outcome_mean,
          random_intercept_sd = scenario$random_intercept_sd_value,
          random_slope_sd = scenario$random_slope_sd_value,
          residual_sd = scenario$residual_sd_value,
          trial_missingness_range = trial_missingness_range,
          subject_dropout_range = subject_dropout_range
        )

        round_n <- sum(vapply(chunk_results, function(x) as.integer(x$chunk_nsim), integer(1)))
        round_rejects <- sum(vapply(chunk_results, function(x) as.integer(x$rejects), integer(1)))
        total_nsim <- total_nsim + round_n
        total_rejects <- total_rejects + round_rejects
      }

      power <- total_rejects / total_nsim
      se <- sqrt(power * (1 - power) / total_nsim)
      out <- data.frame(
        scenario_label = scenario$scenario_label,
        varied_component = scenario$varied_component,
        multiplier = scenario$multiplier,
        random_intercept_sd = scenario$random_intercept_sd_value,
        random_slope_sd = scenario$random_slope_sd_value,
        residual_sd = scenario$residual_sd_value,
        n_subjects = n_subjects,
        power = power,
        lower = pmax(0, power - 1.96 * se),
        upper = pmin(1, power + 1.96 * se),
        nsim = total_nsim,
        stringsAsFactors = FALSE
      )
      cat(sprintf("[%s] DONE scenario=%s N=%d power=%.3f\n", format(Sys.time(), "%H:%M:%S"), scenario$scenario_label, n_subjects, power))
      flush.console()
      out
    })
    do.call(rbind, n_rows)
  })

  power_df <- do.call(rbind, scenario_rows)
  power_df$scenario_label <- factor(power_df$scenario_label, levels = scenario_order, ordered = TRUE)
  power_df$meets_target_90 <- power_df$power >= strict_power_target
  write.csv(power_df, file.path(output_dir, paste0(output_prefix, "_curve.csv")), row.names = FALSE)

  # Summarize minimum sample size for target power
  summary_rows <- do.call(rbind, lapply(split(power_df, power_df$scenario_label), function(df) {
    df <- df[order(df$n_subjects), , drop = FALSE]
    hit <- df[df$meets_target_90, "n_subjects", drop = TRUE]
    data.frame(
      scenario_label = as.character(df$scenario_label[1]),
      varied_component = as.character(df$varied_component[1]),
      multiplier = as.numeric(df$multiplier[1]),
      N_min_for_90 = if (length(hit) > 0) min(hit) else NA_integer_,
      power_at_max_N = df$power[which.max(df$n_subjects)],
      monotonic_non_decreasing = all(diff(df$power) >= -0.02),
      stringsAsFactors = FALSE
    )
  }))
  write.csv(summary_rows, file.path(output_dir, paste0(output_prefix, "_summary.csv")), row.names = FALSE)

  curve_plot <- ggplot(power_df, aes(x = .data$n_subjects, y = .data$power, color = .data$scenario_label, group = .data$scenario_label)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.5) +
    geom_hline(yintercept = strict_power_target, linetype = "dashed", color = "grey40", linewidth = 0.5) +
    scale_x_continuous(breaks = sort(unique(subject_breaks))) +
    scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1), breaks = c(0, 0.25, 0.50, 0.75, 0.90, 1)) +
    labs(x = "Subjects", y = "Power", color = "Scenario", title = plot_title) +
    theme_minimal(base_size = 13) +
    theme(panel.grid.minor = element_blank())
  png(file = file.path(output_dir, paste0(output_prefix, ".png")), width = 2200, height = 1400, res = 300)
  print(curve_plot)
  dev.off()

  # Plot heatmap of power by sample size and scenario
  heatmap_plot <- ggplot(power_df, aes(x = factor(.data$n_subjects), y = .data$scenario_label, fill = .data$power)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.2f", .data$power)), color = "white", size = 3) +
    scale_fill_gradient(low = heat_low_color, high = heat_high_color, limits = c(0, 1)) +
    labs(x = "Subjects", y = "Scenario", fill = "Power", title = paste0(plot_title, " Heatmap")) +
    theme_minimal(base_size = 13) +
    theme(panel.grid = element_blank())
  png(file = file.path(output_dir, paste0(output_prefix, "_heatmap.png")), width = 2200, height = 1400, res = 300)
  print(heatmap_plot)
  dev.off()

  power_df
}

result <- runSESOI()
print(result)
