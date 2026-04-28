## Simulates SESOI-based power for gamma frequency with mixed models and saves a CSV power curve plus PNG figure

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
  seed <- 123L
  alpha <- 0.05
  nsim <- 5000L
  parallel_workers <- 8L
  parallel_round_chunk_nsim <- 1L
  strict_power_target <- 0.90
  subject_breaks <- seq(5L, 75L, by = 5L)
  contrast_levels <- c("25", "50", "75", "100")
  sim_col <- "gamma_frequency"
  target_term <- "contrast_num_c"
  model_formula <- gamma_frequency ~ contrast_num_c + (1 + contrast_num_c | Subject)

  gcp_root <- resolve_gcp_root()
  input_file <- file.path(gcp_root, "data", "features", "GCP_eeg_GED_gamma_metrics_trials.csv")
  output_dir <- file.path(gcp_root, "figures", "power_analysis")
  output_prefix <- "GCP_power_analysis_gamma_frequency_SESOI_only"
  plot_title <- "Power Analysis: Gamma Frequency"

  set.seed(seed)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  observed_trials <- read.csv(input_file, stringsAsFactors = FALSE)
  trial_counts <- aggregate(rep(1, nrow(observed_trials)) ~ Subject + Contrast, data = observed_trials, FUN = length)
  names(trial_counts)[3] <- "n_trials"
  trials_per_condition <- as.integer(round(stats::median(trial_counts$n_trials, na.rm = TRUE)))

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
      Subject <- factor(rep(seq_len(n_subjects), each = length(contrast_levels) * trials_per_condition))
      contrast <- factor(
        rep(rep(contrast_levels, each = trials_per_condition), times = n_subjects),
        levels = contrast_levels,
        ordered = TRUE
      )
      contrast_num <- as.numeric(as.character(contrast))
      contrast_num_c <- as.numeric(scale(contrast_num, center = TRUE, scale = TRUE))
      dat <- data.frame(Subject = Subject, contrast_num_c = contrast_num_c)

      n_subject_levels <- nlevels(dat$Subject)
      ri_sd <- as.numeric(scenario$random_intercept_sd[1])
      rs_sd <- as.numeric(scenario$random_slope_sd[1])
      e_sd <- as.numeric(scenario$residual_sd[1])
      beta_raw <- as.numeric(scenario$beta_raw[1])
      mu0 <- as.numeric(scenario$outcome_mean[1])
      ri_mult <- as.numeric(scenario$random_intercept_sd_multiplier[1])
      rs_mult <- as.numeric(scenario$random_slope_sd_multiplier[1])
      e_mult <- as.numeric(scenario$residual_sd_multiplier[1])

      if (!is.finite(ri_mult)) ri_mult <- 1
      if (!is.finite(rs_mult)) rs_mult <- 1
      if (!is.finite(e_mult)) e_mult <- 1

      random_intercepts <- rnorm(n_subject_levels, mean = 0, sd = ri_sd * ri_mult)
      random_slopes <- rnorm(n_subject_levels, mean = 0, sd = rs_sd * rs_mult)
      x <- dat$contrast_num_c
      mu <- mu0 + random_intercepts[dat$Subject] + random_slopes[dat$Subject] * x + beta_raw * x
      dat[[sim_col]] <- mu + rnorm(nrow(dat), mean = 0, sd = e_sd * e_mult)

      fit <- suppressMessages(lmer(model_formula, data = dat, REML = FALSE))
      cf <- as.data.frame(summary(fit)$coefficients)
      cf$term <- rownames(cf)
      target_row <- cf[cf$term == target_term, , drop = FALSE]
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

  rows <- lapply(subject_breaks, function(n_subjects) {
    cat(sprintf("[%s] START N=%d\n", format(Sys.time(), "%H:%M:%S"), n_subjects))
    flush.console()

    total_nsim <- 0L
    total_rejects <- 0L
    round_idx <- 0L

    while (total_nsim < nsim) {
      remaining <- nsim - total_nsim
      workers_this_round <- min(length(cl), remaining)
      chunk_sizes <- rep(parallel_round_chunk_nsim, workers_this_round)
      overflow <- (workers_this_round * parallel_round_chunk_nsim) - remaining
      if (overflow > 0L) {
        for (k in seq_len(overflow)) {
          idx <- ((k - 1L) %% workers_this_round) + 1L
          chunk_sizes[idx] <- chunk_sizes[idx] - 1L
        }
      }
      chunk_sizes <- chunk_sizes[chunk_sizes > 0L]
      round_idx <- round_idx + 1L

      chunk_results <- parallel::parLapply(
        cl = cl,
        X = chunk_sizes,
        fun = estimate_power_chunk,
        scenario = scenario,
        n_subjects = n_subjects,
        trials_per_condition = trials_per_condition,
        sim_col = sim_col,
        model_formula = model_formula,
        target_term = target_term,
        alpha = alpha
      )

      round_n <- sum(vapply(chunk_results, function(x) as.integer(x$chunk_nsim), integer(1)))
      round_rejects <- sum(vapply(chunk_results, function(x) as.integer(x$rejects), integer(1)))
      total_nsim <- total_nsim + round_n
      total_rejects <- total_rejects + round_rejects
      interim_power <- total_rejects / total_nsim

      cat(sprintf("[%s] N=%d | round=%d | processed=%d/%d | interim_power=%.3f\n", format(Sys.time(), "%H:%M:%S"), n_subjects, round_idx, total_nsim, nsim, interim_power))
      flush.console()
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
    cat(sprintf("[%s] DONE N=%d | power=%.3f\n", format(Sys.time(), "%H:%M:%S"), n_subjects, power))
    flush.console()
    out
  })

  power_df <- do.call(rbind, rows)
  power_df$meets_target_90 <- power_df$power >= strict_power_target
  write.csv(power_df, file.path(output_dir, paste0(output_prefix, "_curve.csv")), row.names = FALSE)

  x_breaks <- sort(unique(subject_breaks))
  plot_obj <- ggplot(power_df, aes(x = .data$n_subjects, y = .data$power)) +
    geom_point(size = 2, color = "blue") +
    geom_line(color = "blue", linetype = "dotted") +
    geom_errorbar(aes(ymin = .data$lower, ymax = .data$upper), width = 3, color = "blue", alpha = 0.5) +
    geom_hline(yintercept = strict_power_target, linetype = "dashed", color = "grey", linewidth = 0.5) +
    scale_x_continuous(breaks = x_breaks) +
    scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1), breaks = c(0, 0.25, 0.50, 0.75, 0.90, 1)) +
    labs(x = "Subjects", y = "Power", title = plot_title) +
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

  png(file = file.path(output_dir, paste0(output_prefix, ".png")), width = 2200, height = 1400, res = 300)
  print(plot_obj)
  dev.off()

  power_df
}

result <- runSESOI()
print(result)
## Simulates SESOI-based power for gamma frequency with mixed models

# install.packages("lme4", repos = "https://cloud.r-project.org") # Install mixed-model package explicitly
# install.packages("ggplot2", repos = "https://cloud.r-project.org") # Install plotting package explicitly
# install.packages("scales", repos = "https://cloud.r-project.org") # Install axis-label helper package explicitly

suppressPackageStartupMessages(library(lme4)) # Load mixed-model engine used in each simulation
suppressPackageStartupMessages(library(ggplot2)) # Load grammar-of-graphics plotting tools
suppressPackageStartupMessages(library(scales)) # Load percent formatter for y-axis labels

# Execute all simulation and plotting steps from one clearly named entry point
runSESOI <- function() {
  # ----------------------------- Fixed study settings -----------------------------
  seed <- 123L # Set deterministic seed to make simulation output reproducible
  alpha <- 0.05 # Set significance threshold used for binary reject/no-reject decision
  nsim <- 5000L # Set number of Monte-Carlo repetitions per sample-size condition
  parallel_workers <- 8L # Set number of parallel PSOCK workers for faster simulation
  parallel_round_chunk_nsim <- 1L # Set minimal chunk size to report progress frequently
  strict_power_target <- 0.90 # Define power criterion used in the output summary column
  subject_breaks <- seq(5L, 75L, by = 5L) # Define tested total subject counts on x-axis
  contrast_levels <- c("25", "50", "75", "100") # Define contrast conditions per subject
  sim_col <- "gamma_frequency" # Define simulated dependent-variable column name
  target_term <- "contrast_num_c" # Define fixed effect term tested for significance
  model_formula <- gamma_frequency ~ contrast_num_c + (1 + contrast_num_c | Subject) # Define fitted LMM for each synthetic dataset

  # ----------------------------- OS-aware project paths -----------------------------
  gcp_root <- resolve_gcp_root() # Resolve shared project root for local macOS and Windows server
  input_file <- file.path(gcp_root, "data", "features", "GCP_eeg_GED_gamma_metrics_trials.csv") # Use one explicit source file for trial counts
  output_dir <- file.path(gcp_root, "figures", "power_analysis") # Use one explicit output folder for CSV and PNG artifacts
  output_prefix <- "GCP_power_analysis_gamma_frequency_SESOI_only" # Use one explicit filename stem for all result files
  plot_title <- "Power Curve: Gamma Frequency SESOI (contrast_num_c)" # Use one explicit title for figure output

  # ----------------------------- Precompute shared quantities -----------------------------
  set.seed(seed) # Initialize random number generator on master process
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE) # Ensure output directory exists before writing files

  observed_trials <- read.csv(input_file, stringsAsFactors = FALSE) # Load empirical per-trial table to infer trials per condition
  trial_counts <- aggregate(rep(1, nrow(observed_trials)) ~ Subject + Contrast, data = observed_trials, FUN = length) # Count empirical trial repetitions
  names(trial_counts)[3] <- "n_trials" # Rename aggregate count column for clarity in downstream median call
  trials_per_condition <- as.integer(round(stats::median(trial_counts$n_trials, na.rm = TRUE))) # Set simulation trial count from observed median

  scenario <- data.frame( # Define SESOI-generating parameters as a one-row data frame
    scenario_label = "SESOI", # Attach readable scenario label to output rows
    scenario_role = "sesoi", # Attach machine-oriented scenario role name
    beta_raw = 0.05, # Set SESOI effect size on scaled contrast predictor
    outcome_mean = 0.00, # Set global mean outcome intercept before random components
    random_intercept_sd = 0.20, # Set between-subject intercept variability
    random_slope_sd = 0.10, # Set between-subject slope variability
    residual_sd = 1.00, # Set within-observation residual noise variability
    random_intercept_sd_multiplier = 1.00, # Keep random-intercept multiplier at baseline
    random_slope_sd_multiplier = 1.00, # Keep random-slope multiplier at baseline
    residual_sd_multiplier = 1.00, # Keep residual multiplier at baseline
    stringsAsFactors = FALSE # Prevent automatic factor conversion for robust numeric extraction
  )

  # ----------------------------- Worker function setup -----------------------------
  # Keep one small worker function because PSOCK workers require exported callable code
  estimate_power_chunk <- function(
      chunk_nsim,
      scenario,
      n_subjects,
      trials_per_condition,
      sim_col,
      model_formula,
      target_term,
      alpha) {
    rejects <- 0L # Initialize reject counter for this chunk
    for (i in seq_len(chunk_nsim)) { # Repeat simulation-fit-test cycle chunk_nsim times
      Subject <- factor(rep(seq_len(n_subjects), each = length(contrast_levels) * trials_per_condition)) # Build grouped subject identifier vector
      contrast <- factor( # Build ordered contrast factor per subject across all trials
        rep(rep(contrast_levels, each = trials_per_condition), times = n_subjects),
        levels = contrast_levels,
        ordered = TRUE
      )
      contrast_num <- as.numeric(as.character(contrast)) # Convert ordered contrast labels to numeric values
      contrast_num_c <- as.numeric(scale(contrast_num, center = TRUE, scale = TRUE)) # Z-scale predictor to stabilize coefficient interpretation
      dat <- data.frame(Subject = Subject, contrast_num_c = contrast_num_c) # Build model-data frame structure

      n_subject_levels <- nlevels(dat$Subject) # Cache number of grouped levels for random effects
      ri_sd <- as.numeric(scenario$random_intercept_sd[1]) # Extract random-intercept SD
      rs_sd <- as.numeric(scenario$random_slope_sd[1]) # Extract random-slope SD
      e_sd <- as.numeric(scenario$residual_sd[1]) # Extract residual SD
      beta_raw <- as.numeric(scenario$beta_raw[1]) # Extract true fixed-effect slope
      mu0 <- as.numeric(scenario$outcome_mean[1]) # Extract global baseline mean
      ri_mult <- as.numeric(scenario$random_intercept_sd_multiplier[1]) # Extract multiplier for intercept SD
      rs_mult <- as.numeric(scenario$random_slope_sd_multiplier[1]) # Extract multiplier for slope SD
      e_mult <- as.numeric(scenario$residual_sd_multiplier[1]) # Extract multiplier for residual SD

      if (!is.finite(ri_mult)) ri_mult <- 1 # Protect against invalid multiplier values
      if (!is.finite(rs_mult)) rs_mult <- 1 # Protect against invalid multiplier values
      if (!is.finite(e_mult)) e_mult <- 1 # Protect against invalid multiplier values

      random_intercepts <- rnorm(n_subject_levels, mean = 0, sd = ri_sd * ri_mult) # Draw subject-specific intercept deviations
      random_slopes <- rnorm(n_subject_levels, mean = 0, sd = rs_sd * rs_mult) # Draw subject-specific slope deviations
      x <- dat$contrast_num_c # Reference scaled predictor once for readability
      mu <- mu0 + random_intercepts[dat$Subject] + random_slopes[dat$Subject] * x + beta_raw * x # Construct expected value per observation
      dat[[sim_col]] <- mu + rnorm(nrow(dat), mean = 0, sd = e_sd * e_mult) # Add residual noise and store simulated outcome

      fit <- suppressMessages(lmer(model_formula, data = dat, REML = FALSE)) # Fit target mixed model on simulated data
      cf <- as.data.frame(summary(fit)$coefficients) # Extract fixed-effect summary table
      cf$term <- rownames(cf) # Promote coefficient names to explicit column for robust indexing
      target_row <- cf[cf$term == target_term, , drop = FALSE] # Isolate target term row for p-value logic
      if ("Pr(>|t|)" %in% names(target_row)) { # Use reported p-value when present
        p_value <- as.numeric(target_row[["Pr(>|t|)"]]) # Read p-value directly from model summary
      } else {
        p_value <- 2 * stats::pnorm(abs(as.numeric(target_row[["t value"]])), lower.tail = FALSE) # Approximate two-sided p-value from z/t statistic
      }
      rejects <- rejects + as.integer(is.finite(p_value) && p_value < alpha) # Update reject count with binary decision
    }
    list(chunk_nsim = as.integer(chunk_nsim), rejects = as.integer(rejects)) # Return chunk statistics for master aggregation
  }

  # ----------------------------- Parallel cluster setup -----------------------------
  cl <- parallel::makeCluster(as.integer(parallel_workers), type = "PSOCK") # Create PSOCK cluster with requested worker count
  on.exit(parallel::stopCluster(cl), add = TRUE) # Guarantee cluster shutdown even if execution errors occur
  parallel::clusterEvalQ(cl, suppressPackageStartupMessages(library(lme4))) # Load lme4 on each worker process
  parallel::clusterSetRNGStream(cl, iseed = seed) # Set reproducible and independent RNG streams across workers
  parallel::clusterExport(cl, varlist = c("estimate_power_chunk", "contrast_levels"), envir = environment()) # Export local worker function and shared contrast settings

  # ----------------------------- Main power-curve loop -----------------------------
  rows <- lapply(subject_breaks, function(n_subjects) { # Iterate over each candidate sample size
    cat(sprintf("[%s] START N=%d\n", format(Sys.time(), "%H:%M:%S"), n_subjects)) # Emit start marker for this sample-size block
    flush.console() # Force progress output in interactive sessions

    total_nsim <- 0L # Track cumulative number of simulations completed for this N
    total_rejects <- 0L # Track cumulative number of significant tests for this N
    round_idx <- 0L # Track round number for progress reporting

    while (total_nsim < nsim) { # Continue until requested simulation budget is exhausted
      remaining <- nsim - total_nsim # Calculate simulations still needed in this N block
      workers_this_round <- min(length(cl), remaining) # Avoid assigning empty work to workers
      chunk_sizes <- rep(parallel_round_chunk_nsim, workers_this_round) # Allocate equal chunk sizes initially
      overflow <- (workers_this_round * parallel_round_chunk_nsim) - remaining # Compute excess jobs beyond remaining budget
      if (overflow > 0L) { # Rebalance chunks when final round would overshoot nsim
        for (k in seq_len(overflow)) { # Decrement chunk sizes cyclically to remove overflow
          idx <- ((k - 1L) %% workers_this_round) + 1L # Rotate subtraction index over active workers
          chunk_sizes[idx] <- chunk_sizes[idx] - 1L # Remove one simulation unit from selected chunk
        }
      }
      chunk_sizes <- chunk_sizes[chunk_sizes > 0L] # Remove any zero-sized chunks after balancing
      round_idx <- round_idx + 1L # Increment round counter before dispatch

      chunk_results <- parallel::parLapply( # Dispatch chunk jobs in parallel to cluster workers
        cl = cl,
        X = chunk_sizes,
        fun = estimate_power_chunk,
        scenario = scenario,
        n_subjects = n_subjects,
        trials_per_condition = trials_per_condition,
        sim_col = sim_col,
        model_formula = model_formula,
        target_term = target_term,
        alpha = alpha
      )

      round_n <- sum(vapply(chunk_results, function(x) as.integer(x$chunk_nsim), integer(1))) # Count completed simulations this round
      round_rejects <- sum(vapply(chunk_results, function(x) as.integer(x$rejects), integer(1))) # Count significant tests this round
      total_nsim <- total_nsim + round_n # Accumulate simulation count
      total_rejects <- total_rejects + round_rejects # Accumulate rejection count
      interim_power <- total_rejects / total_nsim # Compute running power estimate for progress display

      cat(sprintf("[%s] N=%d | round=%d | processed=%d/%d | interim_power=%.3f\n", format(Sys.time(), "%H:%M:%S"), n_subjects, round_idx, total_nsim, nsim, interim_power)) # Report simulation progress
      flush.console() # Flush output to keep console responsive during long runs
    }

    power <- total_rejects / total_nsim # Compute final power estimate at this sample size
    se <- sqrt(power * (1 - power) / total_nsim) # Compute Monte-Carlo standard error for confidence interval
    out <- data.frame( # Build one-row result table for this sample size
      scenario_label = scenario$scenario_label, # Preserve scenario label in output
      scenario_role = scenario$scenario_role, # Preserve scenario role in output
      n_subjects = n_subjects, # Store tested sample size
      power = power, # Store estimated power
      lower = pmax(0, power - 1.96 * se), # Store lower 95% confidence bound
      upper = pmin(1, power + 1.96 * se), # Store upper 95% confidence bound
      nsim = total_nsim # Store realized simulation count
    )
    cat(sprintf("[%s] DONE N=%d | power=%.3f\n", format(Sys.time(), "%H:%M:%S"), n_subjects, power)) # Emit completion marker for this N
    flush.console() # Flush completion marker immediately
    out # Return one-row output for row binding
  })

  # ----------------------------- Save outputs -----------------------------
  power_df <- do.call(rbind, rows) # Stack all per-N rows into a single data frame
  power_df$meets_target_90 <- power_df$power >= strict_power_target # Add logical indicator for target-power attainment
  write.csv(power_df, file.path(output_dir, paste0(output_prefix, "_curve.csv")), row.names = FALSE) # Save tabular power curve

  x_breaks <- sort(unique(subject_breaks)) # Generate explicit x-axis ticks from tested sample sizes
  plot_obj <- ggplot(power_df, aes(x = .data$n_subjects, y = .data$power)) + # Initialize ggplot mapping for power curve
    geom_point(size = 2, color = "blue") + # Draw observed power points
    geom_line(color = "blue", linetype = "dotted") + # Connect points with dotted line for trend visibility
    geom_errorbar(aes(ymin = .data$lower, ymax = .data$upper), width = 3, color = "blue", alpha = 0.5) + # Add confidence interval bars
    geom_hline(yintercept = strict_power_target, linetype = "dashed", color = "grey", linewidth = 0.5) + # Draw horizontal target-power reference line
    scale_x_continuous(breaks = x_breaks) + # Fix x-axis ticks to simulated sample sizes
    scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1), breaks = c(0, 0.25, 0.50, 0.75, 0.90, 1)) + # Show y-axis in percentage format with fixed bounds
    labs(x = "Subjects", y = "Power", title = plot_title) + # Set axis labels and figure title
    theme_minimal(base_size = 15) + # Apply minimal theme with readable baseline font size
    theme( # Override theme details for publication-style frame and ticks
      plot.background = element_rect(fill = "white", colour = NA), # Keep outer plot background white
      panel.background = element_rect(fill = "white", colour = NA), # Keep panel background white
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1), # Draw black panel border
      axis.line = element_line(colour = "black", linewidth = 1), # Draw black axis lines
      panel.grid.major = element_blank(), # Remove major grid lines
      panel.grid.minor = element_blank(), # Remove minor grid lines
      axis.ticks.length = grid::unit(0.2, "cm"), # Set visible tick length
      axis.ticks = element_line(linewidth = 0.5), # Set base tick thickness
      axis.ticks.x = element_line(colour = "black", linewidth = 0.5, lineend = "square"), # Set x tick style
      axis.ticks.y = element_line(colour = "black", linewidth = 0.5, lineend = "square"), # Set y tick style
      axis.text.x = element_text(margin = margin(t = 10)), # Separate x labels from axis line
      axis.text.y = element_text(margin = margin(r = 10)) # Separate y labels from axis line
    )

  png(file = file.path(output_dir, paste0(output_prefix, ".png")), width = 2200, height = 1400, res = 300) # Open PNG device for high-resolution figure export
  print(plot_obj) # Render plot object onto active PNG device
  dev.off() # Close graphics device to finalize image file on disk

  power_df # Return final data frame for console inspection and downstream reuse
}

result <- runSESOI() # Execute SESOI simulation workflow
print(result) # Print final table to console for immediate review
