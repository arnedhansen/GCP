## GCP Pilot Statistics for Power Parameterization
## This script exports pilot descriptive statistics for simulation-based power analyses.

ensure_packages <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    install.packages(missing, repos = "https://cloud.r-project.org")
  }
}

ensure_packages(c("lme4"))
suppressPackageStartupMessages({
  library(lme4)
})

set.seed(123)

resolve_project_roots <- function() {
  if (.Platform$OS.type == "windows") {
    list(
      data_root = file.path("W:/Students/Arne/GCP"),
      share_root = "W:/"
    )
  } else {
    list(
      data_root = file.path("/Volumes/g_psyplafor_methlab$/Students/Arne/GCP"),
      share_root = "/Volumes/g_psyplafor_methlab$/"
    )
  }
}

roots <- resolve_project_roots()
PREFERRED_OUTPUT_DIR <- file.path(roots$data_root, "figures", "power_analysis", "pilot_stats")
FALLBACK_OUTPUT_DIR <- file.path(getwd(), "_power_analysis", "outputs", "pilot_stats")

resolve_input_file <- function() {
  candidates <- c(
    file.path(roots$data_root, "data", "features", "GCP_eeg_GED_gamma_metrics_trials.csv"),
    file.path(getwd(), "data", "features", "GCP_eeg_GED_gamma_metrics_trials.csv")
  )
  hit <- candidates[file.exists(candidates)]
  if (length(hit) == 0) {
    stop(
      "Could not find GCP_eeg_GED_gamma_metrics_trials.csv. Checked: ",
      paste(candidates, collapse = " | ")
    )
  }
  hit[1]
}
INPUT_FILE <- resolve_input_file()

resolve_output_dir <- function() {
  if (dir.exists(dirname(PREFERRED_OUTPUT_DIR))) {
    return(PREFERRED_OUTPUT_DIR)
  }
  FALLBACK_OUTPUT_DIR
}

summary_stats <- function(x) {
  x <- x[is.finite(x)]
  n <- length(x)
  if (n == 0) {
    return(data.frame(
      n = 0, mean = NA_real_, sd = NA_real_, median = NA_real_,
      iqr = NA_real_, se = NA_real_, ci95_low = NA_real_, ci95_high = NA_real_
    ))
  }
  mu <- mean(x)
  s <- stats::sd(x)
  se <- s / sqrt(n)
  ci_low <- mu - 1.96 * se
  ci_high <- mu + 1.96 * se
  data.frame(
    n = n,
    mean = mu,
    sd = s,
    median = stats::median(x),
    iqr = stats::IQR(x),
    se = se,
    ci95_low = ci_low,
    ci95_high = ci_high
  )
}

paired_contrasts <- function(df_subject_means, value_col) {
  cond_levels <- levels(df_subject_means$Contrast)
  pairs <- utils::combn(cond_levels, 2, simplify = FALSE)
  out <- vector("list", length(pairs))

  for (i in seq_along(pairs)) {
    c1 <- pairs[[i]][1]
    c2 <- pairs[[i]][2]

    d1 <- df_subject_means[df_subject_means$Contrast == c1, c("Subject", value_col)]
    d2 <- df_subject_means[df_subject_means$Contrast == c2, c("Subject", value_col)]
    names(d1)[2] <- "v1"
    names(d2)[2] <- "v2"
    m <- merge(d1, d2, by = "Subject", all = FALSE)
    diff_vec <- m$v2 - m$v1

    test_obj <- stats::t.test(m$v2, m$v1, paired = TRUE)
    dz <- mean(diff_vec, na.rm = TRUE) / stats::sd(diff_vec, na.rm = TRUE)

    out[[i]] <- data.frame(
      outcome = value_col,
      contrast_1 = c1,
      contrast_2 = c2,
      n_subjects = nrow(m),
      mean_diff = mean(diff_vec, na.rm = TRUE),
      sd_diff = stats::sd(diff_vec, na.rm = TRUE),
      t_value = unname(test_obj$statistic),
      df = unname(test_obj$parameter),
      p_value = unname(test_obj$p.value),
      ci95_low = unname(test_obj$conf.int[1]),
      ci95_high = unname(test_obj$conf.int[2]),
      cohens_dz = dz
    )
  }

  out_df <- do.call(rbind, out)
  out_df$p_value_bh <- stats::p.adjust(out_df$p_value, method = "BH")
  out_df
}

extract_mixed_model_stats <- function(model_obj, outcome_label) {
  sm <- summary(model_obj)
  fix <- as.data.frame(sm$coefficients)
  fix$term <- rownames(fix)
  rownames(fix) <- NULL

  estimate_col <- grep("^Estimate$", names(fix), value = TRUE, ignore.case = TRUE)[1]
  se_col <- grep("Std\\.?\\s*Error", names(fix), value = TRUE, ignore.case = TRUE)[1]
  t_col <- grep("^t\\s*value$", names(fix), value = TRUE, ignore.case = TRUE)[1]
  z_col <- grep("^z\\s*value$", names(fix), value = TRUE, ignore.case = TRUE)[1]
  df_col <- grep("^df$", names(fix), value = TRUE, ignore.case = TRUE)[1]
  p_col <- grep("^Pr\\(", names(fix), value = TRUE)[1]

  if (is.na(estimate_col) || is.na(se_col)) {
    stop("Could not identify estimate/std_error columns in model summary.")
  }

  stat_col <- if (!is.na(t_col)) t_col else z_col
  if (is.na(stat_col)) {
    stop("Could not identify t/z statistic column in model summary.")
  }

  out_cols <- data.frame(
    term = fix$term,
    estimate = as.numeric(fix[[estimate_col]]),
    std_error = as.numeric(fix[[se_col]]),
    statistic = as.numeric(fix[[stat_col]]),
    stringsAsFactors = FALSE
  )
  if (!is.na(df_col)) {
    out_cols$df <- as.numeric(fix[[df_col]])
  }
  if (!is.na(p_col)) {
    out_cols$p_value <- as.numeric(fix[[p_col]])
  }
  out_cols$outcome <- outcome_label

  ci <- suppressMessages(confint(model_obj, parm = "beta_", method = "Wald"))
  ci_df <- data.frame(
    term = rownames(ci),
    ci95_low = ci[, 1],
    ci95_high = ci[, 2],
    stringsAsFactors = FALSE
  )
  ci_df$term <- sub("^beta_", "", ci_df$term)

  merge(out_cols, ci_df, by = "term", all.x = TRUE)
}

safe_random_slope_sd <- function(model_obj, term_name) {
  vc <- as.data.frame(VarCorr(model_obj))
  hit <- vc[vc$grp == "Subject" & vc$var1 == term_name & vc$var2 == "", "sdcor"]
  if (length(hit) == 0) {
    return(0)
  }
  as.numeric(hit[1])
}

build_triangulated_manifest <- function(subject_level_means, mm_freq, mm_power) {
  ## External assumptions should be reviewed and replaced with published values.
  ## They are intentionally conservative defaults, not pilot-derived anchors.
  external_effects <- list(
    PeakFrequency = list(
      sesoi_std = 0.10,
      external_lower_std = 0.12,
      external_point_std = 0.18
    ),
    PeakAmplitude = list(
      sesoi_std = -0.08,
      external_lower_std = -0.10,
      external_point_std = -0.16
    )
  )

  make_rows <- function(outcome, target_term, model_obj) {
    y <- subject_level_means[[outcome]]
    y_sd <- stats::sd(y, na.rm = TRUE)
    y_mean <- mean(y, na.rm = TRUE)
    pilot_fix <- lme4::fixef(model_obj)
    pilot_beta <- unname(pilot_fix[target_term])
    pilot_beta_std <- pilot_beta / y_sd
    ext <- external_effects[[outcome]]

    random_intercept_sd <- as.numeric(attr(VarCorr(model_obj)$Subject, "stddev")[1])
    random_slope_sd <- safe_random_slope_sd(model_obj, "contrast_num_c")
    residual_sd <- sigma(model_obj)

    data.frame(
      outcome = outcome,
      target_term = target_term,
      scenario_label = c("Pessimistic", "Base", "SESOI", "PilotSecondary"),
      scenario_role = c("pessimistic", "base", "sesoi", "pilot_secondary"),
      effect_source = c("external_lower_bound", "external_point_estimate", "sesoi", "pilot_secondary"),
      beta_std = c(ext$external_lower_std, ext$external_point_std, ext$sesoi_std, pilot_beta_std),
      beta_raw = c(ext$external_lower_std, ext$external_point_std, ext$sesoi_std, pilot_beta_std) * y_sd,
      outcome_mean = y_mean,
      outcome_sd = y_sd,
      random_intercept_sd = random_intercept_sd,
      random_slope_sd = random_slope_sd,
      residual_sd = residual_sd,
      n_subjects_pilot = length(unique(subject_level_means$Subject)),
      stringsAsFactors = FALSE
    )
  }

  rbind(
    make_rows("PeakFrequency", "contrast_num_c", mm_freq),
    make_rows("PeakAmplitude", "contrast_num_c2", mm_power)
  )
}

message("Loading pilot data from: ", INPUT_FILE)
dat <- read.csv(INPUT_FILE, stringsAsFactors = FALSE)

required_cols <- c("Subject", "Contrast", "PeakFrequency", "PeakAmplitude")
missing_cols <- setdiff(required_cols, names(dat))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

dat$Subject <- factor(dat$Subject)
dat$Contrast <- factor(dat$Contrast, levels = c(25, 50, 75, 100), ordered = TRUE)
dat$contrast_num <- as.numeric(as.character(dat$Contrast))
dat$contrast_num_c <- as.numeric(scale(dat$contrast_num, center = TRUE, scale = TRUE))
dat$contrast_num_c2 <- dat$contrast_num_c^2

message("Computing trial-count diagnostics.")
trial_counts_by_subject_condition <- aggregate(
  rep(1, nrow(dat)) ~ Subject + Contrast,
  data = dat,
  FUN = length
)
names(trial_counts_by_subject_condition)[3] <- "n_trials"

message("Computing trial-level descriptive stats by condition.")
trial_level_frequency <- do.call(
  rbind,
  lapply(levels(dat$Contrast), function(cc) {
    x <- dat$PeakFrequency[dat$Contrast == cc]
    cbind(outcome = "PeakFrequency", contrast = cc, summary_stats(x))
  })
)

trial_level_power <- do.call(
  rbind,
  lapply(levels(dat$Contrast), function(cc) {
    x <- dat$PeakAmplitude[dat$Contrast == cc]
    cbind(outcome = "PeakAmplitude", contrast = cc, summary_stats(x))
  })
)

trial_level_descriptives <- rbind(trial_level_frequency, trial_level_power)

message("Computing subject-level means and descriptive stats.")
subject_level_means <- aggregate(
  cbind(PeakFrequency, PeakAmplitude) ~ Subject + Contrast + contrast_num + contrast_num_c + contrast_num_c2,
  data = dat,
  FUN = mean,
  na.rm = TRUE
)

subject_level_frequency <- do.call(
  rbind,
  lapply(levels(subject_level_means$Contrast), function(cc) {
    x <- subject_level_means$PeakFrequency[subject_level_means$Contrast == cc]
    cbind(outcome = "PeakFrequency", contrast = cc, summary_stats(x))
  })
)

subject_level_power <- do.call(
  rbind,
  lapply(levels(subject_level_means$Contrast), function(cc) {
    x <- subject_level_means$PeakAmplitude[subject_level_means$Contrast == cc]
    cbind(outcome = "PeakAmplitude", contrast = cc, summary_stats(x))
  })
)

subject_level_descriptives <- rbind(subject_level_frequency, subject_level_power)

message("Running paired condition contrasts on subject-level means.")
pairwise_frequency <- paired_contrasts(subject_level_means, "PeakFrequency")
pairwise_power <- paired_contrasts(subject_level_means, "PeakAmplitude")
pairwise_contrasts <- rbind(pairwise_frequency, pairwise_power)

message("Fitting mixed models for parameterization guidance.")
mm_freq <- lmer(
  PeakFrequency ~ contrast_num_c + contrast_num_c2 + (1 + contrast_num_c | Subject),
  data = subject_level_means,
  REML = FALSE
)
if (isSingular(mm_freq, tol = 1e-4)) {
  message("Frequency model singular; using random-intercept fallback.")
  mm_freq <- lmer(
    PeakFrequency ~ contrast_num_c + contrast_num_c2 + (1 | Subject),
    data = subject_level_means,
    REML = FALSE
  )
}

mm_power <- lmer(
  PeakAmplitude ~ contrast_num_c + contrast_num_c2 + (1 + contrast_num_c | Subject),
  data = subject_level_means,
  REML = FALSE
)
if (isSingular(mm_power, tol = 1e-4)) {
  message("Power model singular; using random-intercept fallback.")
  mm_power <- lmer(
    PeakAmplitude ~ contrast_num_c + contrast_num_c2 + (1 | Subject),
    data = subject_level_means,
    REML = FALSE
  )
}

mixed_model_stats <- rbind(
  extract_mixed_model_stats(mm_freq, "PeakFrequency"),
  extract_mixed_model_stats(mm_power, "PeakAmplitude")
)

random_effects_summary <- rbind(
  data.frame(outcome = "PeakFrequency", as.data.frame(VarCorr(mm_freq))),
  data.frame(outcome = "PeakAmplitude", as.data.frame(VarCorr(mm_power)))
)

triangulated_manifest <- build_triangulated_manifest(subject_level_means, mm_freq, mm_power)

output_dir <- resolve_output_dir()
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
message("Saving outputs to: ", output_dir)

write.csv(
  trial_counts_by_subject_condition,
  file.path(output_dir, "pilot_trial_counts_by_subject_condition.csv"),
  row.names = FALSE
)
write.csv(
  trial_level_descriptives,
  file.path(output_dir, "pilot_trial_level_descriptives_gamma_freq_power.csv"),
  row.names = FALSE
)
write.csv(
  subject_level_means,
  file.path(output_dir, "pilot_subject_level_means_gamma_freq_power.csv"),
  row.names = FALSE
)
write.csv(
  subject_level_descriptives,
  file.path(output_dir, "pilot_subject_level_descriptives_gamma_freq_power.csv"),
  row.names = FALSE
)
write.csv(
  pairwise_contrasts,
  file.path(output_dir, "pilot_pairwise_contrasts_gamma_freq_power.csv"),
  row.names = FALSE
)
write.csv(
  mixed_model_stats,
  file.path(output_dir, "pilot_mixed_model_fixed_effects_gamma_freq_power.csv"),
  row.names = FALSE
)
write.csv(
  random_effects_summary,
  file.path(output_dir, "pilot_mixed_model_random_effects_gamma_freq_power.csv"),
  row.names = FALSE
)
write.csv(
  triangulated_manifest,
  file.path(output_dir, "power_parameter_manifest.csv"),
  row.names = FALSE
)

saveRDS(
  list(
    trial_counts_by_subject_condition = trial_counts_by_subject_condition,
    trial_level_descriptives = trial_level_descriptives,
    subject_level_means = subject_level_means,
    subject_level_descriptives = subject_level_descriptives,
    pairwise_contrasts = pairwise_contrasts,
    mixed_model_fixed_effects = mixed_model_stats,
    mixed_model_random_effects = random_effects_summary,
    triangulated_manifest = triangulated_manifest
  ),
  file = file.path(output_dir, "pilot_stats_for_power_analysis.rds")
)

message("Pilot statistics export complete.")
