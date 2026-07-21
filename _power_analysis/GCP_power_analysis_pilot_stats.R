## Pilot variance estimates for the linear gamma power and frequency models.

ensure_packages <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0L) {
    install.packages(missing, repos = "https://cloud.r-project.org")
  }
}

ensure_packages("lme4")
suppressPackageStartupMessages(library(lme4))

script_dir <- {
  cmd <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd, value = TRUE)
  if (length(file_arg) > 0L) {
    dirname(sub("^--file=", "", file_arg[1]))
  } else {
    file.path(getwd(), "_power_analysis")
  }
}
source(file.path(script_dir, "GCP_power_analysis_helpers.R"))
# nolint start: object_usage_linter

resolve_gcp_root_for_pilot <- function() {
  candidates <- if (.Platform$OS.type == "windows") {
    c("W:/Students/Arne/GCP", getwd(), dirname(getwd()))
  } else {
    c(
      "/Volumes/g_psyplafor_methlab$/Students/Arne/GCP",
      getwd(),
      dirname(getwd())
    )
  }
  probe <- file.path(candidates, "data", "features", "GCP_merged_data.csv")
  hit <- candidates[file.exists(probe)]
  if (length(hit) == 0L) {
    stop("Could not locate GCP_merged_data.csv. Checked: ", paste(probe, collapse = " | "))
  }
  hit[1]
}

normalize_contrast <- function(x) {
  value <- suppressWarnings(as.numeric(as.character(x)))
  if (all(value %in% 1:4, na.rm = TRUE)) {
    value <- c(25, 50, 75, 100)[as.integer(value)]
  }
  value
}

prepare_gamma_data <- function(raw) {
  id_col <- intersect(c("ID", "Subject"), names(raw))[1]
  frequency_col <- intersect(
    c("Frequency", "PeakFrequency", "GED_PeakFrequency"),
    names(raw)
  )[1]
  power_col <- intersect(
    c("Power", "PeakAmplitude", "GED_PeakAmplitude"),
    names(raw)
  )[1]
  required <- c(id_col, "Condition", frequency_col, power_col)
  if (any(is.na(required)) || !all(required %in% names(raw))) {
    stop("Master matrix must contain subject, condition, gamma frequency, and gamma power.")
  }

  dat <- data.frame(
    Subject = raw[[id_col]],
    contrast_num = normalize_contrast(raw$Condition),
    PeakFrequency = as.numeric(raw[[frequency_col]]),
    PeakAmplitude = as.numeric(raw[[power_col]]),
    stringsAsFactors = FALSE
  )
  if ("Include" %in% names(raw)) {
    include <- raw$Include
    include <- if (is.logical(include)) include else as.integer(include) == 1L
    dat <- dat[include, , drop = FALSE]
  }
  dat <- dat[
    is.finite(dat$contrast_num) &
      is.finite(dat$PeakFrequency) &
      is.finite(dat$PeakAmplitude),
    ,
    drop = FALSE
  ]
  dat <- stats::aggregate(
    cbind(PeakFrequency, PeakAmplitude) ~ Subject + contrast_num,
    data = dat,
    FUN = mean,
    na.rm = TRUE
  )
  dat$Subject <- factor(dat$Subject)
  population_sd <- sqrt(mean((c(25, 50, 75, 100) - 62.5)^2))
  dat$contrast_num_c <- (dat$contrast_num - 62.5) / population_sd
  assert_one_row_per_subject_condition(dat)
  dat
}

summary_stats <- function(x) {
  x <- x[is.finite(x)]
  data.frame(
    n = length(x),
    mean = mean(x),
    sd = stats::sd(x),
    median = stats::median(x),
    iqr = stats::IQR(x),
    stringsAsFactors = FALSE
  )
}

paired_contrasts <- function(dat, outcome) {
  pairs <- utils::combn(c(25, 50, 75, 100), 2, simplify = FALSE)
  do.call(rbind, lapply(pairs, function(pair) {
    left <- dat[dat$contrast_num == pair[1], c("Subject", outcome)]
    right <- dat[dat$contrast_num == pair[2], c("Subject", outcome)]
    names(left)[2] <- "left"
    names(right)[2] <- "right"
    merged <- merge(left, right, by = "Subject")
    difference <- merged$right - merged$left
    test <- stats::t.test(merged$right, merged$left, paired = TRUE)
    data.frame(
      outcome = outcome,
      contrast_1 = pair[1],
      contrast_2 = pair[2],
      n_subjects = nrow(merged),
      mean_difference = mean(difference),
      sd_difference = stats::sd(difference),
      t_value = unname(test$statistic),
      p_value = test$p.value,
      cohens_dz = mean(difference) / stats::sd(difference),
      stringsAsFactors = FALSE
    )
  }))
}

extract_fixed_effects <- function(fit, outcome) {
  coefs <- as.data.frame(summary(fit)$coefficients)
  coefs$term <- rownames(coefs)
  rownames(coefs) <- NULL
  data.frame(
    outcome = outcome,
    model_scope = "linear",
    term = coefs$term,
    estimate = coefs$Estimate,
    std_error = coefs$`Std. Error`,
    statistic = coefs$`t value`,
    stringsAsFactors = FALSE
  )
}

run_pilot_stats <- function(bootstrap_B = 500L) {
  set.seed(123L)
  gcp_root <- resolve_gcp_root_for_pilot()
  input_path <- file.path(gcp_root, "data", "features", "GCP_merged_data.csv")
  output_dir <- file.path(gcp_root, "data", "pilot_stats")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  dat <- prepare_gamma_data(utils::read.csv(input_path, stringsAsFactors = FALSE))
  if (nlevels(dat$Subject) < 3L) {
    stop("At least three pilot participants are required.")
  }

  frequency_fit_info <- fit_lmer_with_fallbacks(
    build_subject_level_formula_list("PeakFrequency", model_scope = "linear"),
    dat
  )
  power_fit_info <- fit_lmer_with_fallbacks(
    build_subject_level_formula_list("PeakAmplitude", model_scope = "linear"),
    dat
  )
  if (is.null(frequency_fit_info) || is.null(power_fit_info)) {
    stop("At least one linear gamma pilot model failed to fit.")
  }
  frequency_fit <- frequency_fit_info$fit
  power_fit <- power_fit_info$fit

  message("Bootstrapping linear gamma frequency variance partition.")
  frequency_bootstrap <- bootstrap_subject_level_variance_partition(
    dat,
    stats::formula(frequency_fit),
    B = bootstrap_B,
    rng_seed = 123L
  )
  message("Bootstrapping linear gamma power variance partition.")
  power_bootstrap <- bootstrap_subject_level_variance_partition(
    dat,
    stats::formula(power_fit),
    B = bootstrap_B,
    rng_seed = 124L
  )

  variance_summary <- rbind(
    summarize_variance_partition_bootstrap(
      frequency_bootstrap,
      "PeakFrequency",
      "linear",
      basename(input_path),
      nlevels(dat$Subject),
      paste(trimws(deparse(stats::formula(frequency_fit))), collapse = " "),
      outcome_mean = mean(dat$PeakFrequency)
    ),
    summarize_variance_partition_bootstrap(
      power_bootstrap,
      "PeakAmplitude",
      "linear",
      basename(input_path),
      nlevels(dat$Subject),
      paste(trimws(deparse(stats::formula(power_fit))), collapse = " "),
      outcome_mean = mean(dat$PeakAmplitude)
    )
  )

  descriptives <- do.call(rbind, lapply(
    c("PeakFrequency", "PeakAmplitude"),
    function(outcome) {
      do.call(rbind, lapply(c(25, 50, 75, 100), function(contrast) {
        cbind(
          outcome = outcome,
          contrast = contrast,
          summary_stats(dat[dat$contrast_num == contrast, outcome])
        )
      }))
    }
  ))
  pairwise <- rbind(
    paired_contrasts(dat, "PeakFrequency"),
    paired_contrasts(dat, "PeakAmplitude")
  )
  fixed_effects <- rbind(
    extract_fixed_effects(frequency_fit, "PeakFrequency"),
    extract_fixed_effects(power_fit, "PeakAmplitude")
  )
  random_effects <- rbind(
    cbind(
      outcome = "PeakFrequency",
      model_scope = "linear",
      as.data.frame(lme4::VarCorr(frequency_fit))
    ),
    cbind(
      outcome = "PeakAmplitude",
      model_scope = "linear",
      as.data.frame(lme4::VarCorr(power_fit))
    )
  )
  literature_effects <- data.frame(
    outcome = c("PeakFrequency", "PeakAmplitude"),
    model_scope = "linear",
    target_term = "contrast_num_c",
    cohens_dz = c(0.45, 0.52),
    source = "Karvat, Ofir, and Landau (2024)",
    comparison = "100% versus 50% contrast",
    stringsAsFactors = FALSE
  )

  utils::write.csv(
    dat,
    file.path(output_dir, "pilot_subject_level_means_gamma_freq_power.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    descriptives,
    file.path(output_dir, "pilot_subject_level_descriptives_gamma_freq_power.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    pairwise,
    file.path(output_dir, "pilot_pairwise_contrasts_gamma_freq_power.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    fixed_effects,
    file.path(output_dir, "pilot_mixed_model_fixed_effects_gamma_freq_power.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    random_effects,
    file.path(output_dir, "pilot_mixed_model_random_effects_gamma_freq_power.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    variance_summary,
    file.path(output_dir, "pilot_subject_level_variance_partition_summary.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    literature_effects,
    file.path(output_dir, "literature_effect_size_manifest.csv"),
    row.names = FALSE
  )
  saveRDS(
    list(
      subject_level_means = dat,
      descriptives = descriptives,
      pairwise_contrasts = pairwise,
      fixed_effects = fixed_effects,
      random_effects = random_effects,
      variance_partition_summary = variance_summary,
      variance_partition_draws = list(
        PeakFrequency = frequency_bootstrap$draws,
        PeakAmplitude = power_bootstrap$draws
      ),
      literature_effects = literature_effects
    ),
    file.path(output_dir, "pilot_stats_for_power_analysis.rds")
  )

  message("Pilot gamma statistics saved to: ", output_dir)
  invisible(list(
    data = dat,
    frequency_fit = frequency_fit,
    power_fit = power_fit,
    variance_summary = variance_summary
  ))
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  bootstrap_B <- 500L
  bootstrap_arg <- grep("^--bootstrap=", args, value = TRUE)
  if (length(bootstrap_arg) > 0L) {
    bootstrap_B <- as.integer(sub("^--bootstrap=", "", bootstrap_arg[1]))
  }
  run_pilot_stats(bootstrap_B = bootstrap_B)
}
# nolint end
