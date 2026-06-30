## SESOI power simulation for linear contrast (subject-level condition means).
## Residual and random-effects sensitivity axes from pilot variance-partition bootstrap.

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(scales))

engine_path <- file.path(getwd(), "_power_analysis", "GCP_power_analysis_SESOI_engine.R")
if (!file.exists(engine_path)) {
  cmd <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd, value = TRUE)
  if (length(file_arg) > 0) {
    engine_path <- file.path(dirname(sub("^--file=", "", file_arg[1])), "GCP_power_analysis_SESOI_engine.R")
  }
}
if (!file.exists(engine_path)) {
  stop("Could not find GCP_power_analysis_SESOI_engine.R")
}
source(engine_path)

runSESOI <- function(
    plot_only = FALSE,
    nsim = 1000,
    subject_breaks = seq(10, 100, by = 10),
    parallel_workers = 8) {
  config <- list(
    label = "SESOI linear PeakFrequency",
    pilot_outcome = "PeakFrequency",
    model_scope = "linear",
    sim_col = "PeakFrequency",
    target_term = "contrast_num_c",
    sesoi_beta = 0.10,
    true_beta = 1.5 * 0.10,
    linear_nuisance_beta = 0,
    subject_dropout_rate = 0.10,
    default_outcome_mean = 53.75,
    output_prefix_base = "GCP_power_analysis_SESOI_linear",
    plot_title_base = "Power Analysis: SESOI Linear Slope (subject-level condition means)"
  )
  run_sesoi_power(
    config = config,
    plot_only = plot_only,
    nsim = nsim,
    subject_breaks = subject_breaks,
    parallel_workers = parallel_workers,
    sensitivity_axes = c("residual", "random_effects")
  )
}
