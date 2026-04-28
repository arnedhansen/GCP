## GCP simulation-based power analysis for gamma power.
## Uses triangulated assumptions (SESOI + external + pilot-secondary).

resolve_script_dir <- function() {
  cmd <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  }
  cwd <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  if (basename(cwd) == "_power_analysis") {
    return(cwd)
  }
  file.path(cwd, "_power_analysis")
}
source(file.path(resolve_script_dir(), "GCP_power_analysis_common.R"))

result <- run_outcome_power_analysis(list(
  manifest_outcome = "PeakAmplitude",
  sim_col = "gamma_power",
  target_term = "contrast_num_c2",
  model_formula = gamma_power ~ contrast_num_c + contrast_num_c2 + (1 + contrast_num_c | Subject),
  linear_nuisance_beta = 0.05,
  alpha = 0.05,
  run_mode = "exploratory",
  nsim_exploratory = 1000,
  nsim_reporting = 5000,
  verbose = TRUE,
  progress_every = 250,
  parallel_enabled = TRUE,
  parallel_workers = 8,
  validation_nsim = 1000,
  validation_n = 40,
  alpha_tolerance = 0.02,
  large_n_min_power = 0.90,
  strict_power_base = 0.90,
  strict_power_pessimistic = 0.80,
  feasible_subject_cap = 80,
  subject_breaks = seq(20, 140, by = 10),
  trials_per_condition = resolve_trials_per_condition(default_value = 160),
  seed = 124,
  file_prefix = "GCP_power_analysis_gamma_power",
  plot_title = "Power Curves: Gamma Power (contrast_num_c2)"
))

print(result$recommendations)
