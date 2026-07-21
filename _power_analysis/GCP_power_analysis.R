## Literature effect size power analysis for gamma outcomes.
## Both outcomes use linear contrast models with a by-subject random slope,
## simulated with SIMR. The sensitivity axis is the random-slope SD.

## Random-slope SD scenarios, expressed as multiples of the fixed linear slope.
## The between-subject SD of the contrast effect is the dominant term governing
## how power grows with N, so it is varied here instead of the residual SD.

#USEAGE
#setwd("~/GitHub/GCP/_power_analysis")
#setwd("..")
#
#source("_power_analysis/GCP_power_analysis_pilot_stats.R")
#run_pilot_stats()
#
#source("_power_analysis/GCP_power_analysis.R")
#run_literature_power_analysis()
#
#source("_power_analysis/GCP_power_analysis_visualizations.R")
#run_power_visualizations()

RANDOM_SLOPE_SD_MULTIPLIERS <- c(low = 0.5, median = 1.0, high = 1.5)

ensure_packages <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0L) {
    install.packages(missing, repos = "https://cloud.r-project.org")
  }
}

ensure_packages(c("lme4", "simr"))
suppressPackageStartupMessages({
  library(lme4)
  library(simr)
})

script_dir <- {
  cmd <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd, value = TRUE)
  if (length(file_arg) > 0L) {
    dirname(sub("^--file=", "", file_arg[1]))
  } else {
    file.path(getwd(), "_power_analysis")
  }
}

resolve_gcp_root <- function() {
  candidates <- if (.Platform$OS.type == "windows") {
    c("W:/Students/Arne/GCP", getwd(), dirname(getwd()))
  } else {
    c(
      "/Volumes/g_psyplafor_methlab$/Students/Arne/GCP",
      getwd(),
      dirname(getwd())
    )
  }
  probe <- file.path(
    candidates,
    "data",
    "pilot_stats",
    "pilot_subject_level_variance_partition_summary.csv"
  )
  hit <- candidates[file.exists(probe)]
  if (length(hit) == 0L) {
    stop(
      "Could not locate pilot_subject_level_variance_partition_summary.csv. Checked: ",
      paste(probe, collapse = " | ")
    )
  }
  hit[1]
}

load_variance_parameters <- function(gcp_root, outcome, model_scope = "linear") {
  path <- file.path(
    gcp_root,
    "data",
    "pilot_stats",
    "pilot_subject_level_variance_partition_summary.csv"
  )
  tab <- utils::read.csv(path, stringsAsFactors = FALSE)
  meta <- tab[
    tab$record_type == "meta" &
      tab$outcome == outcome &
      tab$model_scope == model_scope,
    ,
    drop = FALSE
  ]
  components <- tab[
    tab$record_type == "component" &
      tab$outcome == outcome &
      tab$model_scope == model_scope,
    ,
    drop = FALSE
  ]
  if (nrow(meta) == 0L && outcome == "PeakAmplitude" && model_scope == "linear") {
    warning(
      "Linear PeakAmplitude variance estimates are unavailable. ",
      "Using the existing quadratic pilot variance partition until ",
      "GCP_power_analysis_pilot_stats.R is rerun.",
      call. = FALSE
    )
    meta <- tab[
      tab$record_type == "meta" &
        tab$outcome == outcome &
        tab$model_scope == "quadratic",
      ,
      drop = FALSE
    ]
    components <- tab[
      tab$record_type == "component" &
        tab$outcome == outcome &
        tab$model_scope == "quadratic",
      ,
      drop = FALSE
    ]
  }
  if (nrow(meta) != 1L) {
    stop("Expected one metadata row for ", outcome, " and model scope ", model_scope, ".")
  }

  get_component <- function(component) {
    row <- components[components$component == component, , drop = FALSE]
    if (nrow(row) != 1L) {
      stop("Expected one ", component, " row for ", outcome, ".")
    }
    stats::setNames(
      as.numeric(row[1, c("q025", "q50", "q975")]),
      c("low", "median", "high")
    )
  }

  list(
    outcome_mean = as.numeric(meta$outcome_mean[1]),
    pilot_n = as.integer(meta$pilot_n_subjects[1]),
    residual_sd = get_component("residual_sd"),
    random_intercept_sd = get_component("random_intercept_sd"),
    source_path = path
  )
}

contrast_code <- function(contrast) {
  contrast_values <- c(25, 50, 75, 100)
  population_sd <- sqrt(mean((contrast_values - mean(contrast_values))^2))
  (as.numeric(contrast) - mean(contrast_values)) / population_sd
}

make_simr_design <- function(n_subjects) {
  contrast <- rep(c(25, 50, 75, 100), times = n_subjects)
  data.frame(
    Subject = factor(rep(seq_len(n_subjects), each = 4L)),
    contrast = contrast,
    contrast_num_c = contrast_code(contrast)
  )
}

dz_to_linear_beta <- function(dz, reference_residual_sd) {
  delta_x <- contrast_code(100) - contrast_code(50)
  paired_mean_difference <- as.numeric(dz) * sqrt(2) * reference_residual_sd
  paired_mean_difference / delta_x
}

## Random-slope SD as a multiple of the fixed slope magnitude. Keeping most
## subjects on the same side of zero requires the multiplier to stay near or
## below 1; larger values imply many subjects with null or reversed effects.
slope_sd_from_beta <- function(beta, multiplier) {
  abs(as.numeric(beta)) * as.numeric(multiplier)
}

make_simr_model <- function(
    n_subjects,
    outcome_mean,
    beta,
    random_intercept_sd,
    random_slope_sd,
    residual_sd) {
  dat <- make_simr_design(n_subjects)
  subject_vcov <- matrix(
    c(random_intercept_sd^2, 0, 0, random_slope_sd^2),
    nrow = 2L,
    dimnames = list(
      c("(Intercept)", "contrast_num_c"),
      c("(Intercept)", "contrast_num_c")
    )
  )
  simr::makeLmer(
    outcome ~ contrast_num_c + (1 + contrast_num_c | Subject),
    fixef = c(outcome_mean, beta),
    VarCorr = list(Subject = subject_vcov),
    sigma = residual_sd,
    data = dat
  )
}

extract_power_curve <- function(power_curve, subject_breaks) {
  do.call(rbind, lapply(seq_along(power_curve$ps), function(i) {
    estimate <- power_curve$ps[[i]]
    ci <- suppressWarnings(stats::confint(estimate))
    data.frame(
      n_subjects = subject_breaks[i],
      power = as.numeric(estimate$x / estimate$n),
      lower = as.numeric(ci[1, 1]),
      upper = as.numeric(ci[1, 2]),
      nsim = as.integer(estimate$n),
      stringsAsFactors = FALSE
    )
  }))
}

save_standard_simr_plot <- function(power_curve, path) {
  grDevices::png(path, width = 2200, height = 1400, res = 300)
  print(plot(power_curve))
  grDevices::dev.off()
}

run_literature_simr_curve <- function(
    label,
    outcome,
    dz,
    output_prefix,
    nsim = 1000L,
    subject_breaks = seq(10, 100, by = 10),
    seed = 123L,
    plot_only = FALSE) {
  gcp_root <- resolve_gcp_root()
  figure_dir <- file.path(gcp_root, "figures", "power_analysis", "backup")
  data_dir <- file.path(gcp_root, "data", "power_analysis")
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

  parameters <- load_variance_parameters(gcp_root, outcome, model_scope = "linear")
  reference_sigma <- unname(parameters$residual_sd["median"])
  beta <- dz_to_linear_beta(dz, reference_sigma)
  reference_slope_sd <- slope_sd_from_beta(beta, RANDOM_SLOPE_SD_MULTIPLIERS[["median"]])
  curve_rds <- file.path(data_dir, paste0(output_prefix, "_curve.rds"))
  curve_csv <- file.path(data_dir, paste0(output_prefix, "_curve.csv"))
  figure_path <- file.path(figure_dir, paste0(output_prefix, ".png"))

  message(sprintf(
    "%s | d_z = %.3f | residual SD = %.4f | linear beta = %.4f | random-slope SD = %.4f",
    label, dz, reference_sigma, beta, reference_slope_sd
  ))

  if (plot_only) {
    if (!file.exists(curve_rds)) {
      stop("Cached SIMR power curve not found: ", curve_rds)
    }
    power_curve <- readRDS(curve_rds)
  } else {
    model <- make_simr_model(
      n_subjects = max(subject_breaks),
      outcome_mean = parameters$outcome_mean,
      beta = beta,
      random_intercept_sd = unname(parameters$random_intercept_sd["median"]),
      random_slope_sd = reference_slope_sd,
      residual_sd = reference_sigma
    )
    set.seed(seed)
    power_curve <- suppressWarnings(simr::powerCurve(
      model,
      test = simr::fixed("contrast_num_c", method = "t"),
      along = "Subject",
      breaks = subject_breaks,
      nsim = nsim,
      progress = TRUE
    ))
    saveRDS(power_curve, curve_rds)
    curve_df <- extract_power_curve(power_curve, subject_breaks)
    curve_df$outcome <- label
    curve_df$cohens_dz <- dz
    curve_df$linear_beta <- beta
    curve_df$reference_residual_sd <- reference_sigma
    curve_df$random_slope_sd <- reference_slope_sd
    utils::write.csv(curve_df, curve_csv, row.names = FALSE)
  }

  save_standard_simr_plot(power_curve, figure_path)
  message("Saved standard SIMR curve: ", figure_path)
  invisible(list(
    curve = power_curve,
    parameters = parameters,
    dz = dz,
    beta = beta,
    figure_path = figure_path
  ))
}

run_literature_power_analysis <- function(
    nsim = 1000L,
    subject_breaks = seq(10, 100, by = 10),
    plot_only = FALSE) {
  frequency <- run_literature_simr_curve(
    label = "Gamma peak frequency",
    outcome = "PeakFrequency",
    dz = 0.45,
    output_prefix = "GCP_power_analysis_gamma_frequency_litES",
    nsim = nsim,
    subject_breaks = subject_breaks,
    seed = 123L,
    plot_only = plot_only
  )
  power <- run_literature_simr_curve(
    label = "Gamma peak power",
    outcome = "PeakAmplitude",
    dz = 0.52,
    output_prefix = "GCP_power_analysis_gamma_power_litES",
    nsim = nsim,
    subject_breaks = subject_breaks,
    seed = 124L,
    plot_only = plot_only
  )
  invisible(list(frequency = frequency, power = power))
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  nsim <- 1000L
  nsim_arg <- grep("^--nsim=", args, value = TRUE)
  if (length(nsim_arg) > 0L) {
    nsim <- as.integer(sub("^--nsim=", "", nsim_arg[1]))
  }
  run_literature_power_analysis(
    nsim = nsim,
    plot_only = "--plot-only" %in% args
  )
}
