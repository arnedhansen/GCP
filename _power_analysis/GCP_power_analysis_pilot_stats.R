## GCP Pilot Statistics for SIMR Power Parameterization
## This script summarizes pilot gamma frequency/power by condition and
## exports statistics that can be used to parameterize simulation-based power analyses.

suppressPackageStartupMessages({
  library(lme4)
})

set.seed(123)

INPUT_FILE <- "/Volumes/g_psyplafor_methlab$/Students/Arne/GCP/data/features/GCP_eeg_GED_gamma_metrics_trials.csv"
PREFERRED_OUTPUT_DIR <- "/Volumes/g_psyplafor_methlab$/Students/Arne/GCP/figures/power_analysis/pilot_stats"
FALLBACK_OUTPUT_DIR <- file.path(getwd(), "_power_analysis", "outputs", "pilot_stats")

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
  names(fix) <- c("estimate", "std_error", "statistic", "term")
  fix <- fix[, c("term", "estimate", "std_error", "statistic")]
  fix$outcome <- outcome_label

  ci <- suppressMessages(confint(model_obj, parm = "beta_", method = "Wald"))
  ci_df <- data.frame(
    term = rownames(ci),
    ci95_low = ci[, 1],
    ci95_high = ci[, 2],
    stringsAsFactors = FALSE
  )
  ci_df$term <- sub("^beta_", "", ci_df$term)

  merge(fix, ci_df, by = "term", all.x = TRUE)
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
  Trial ~ Subject + Contrast,
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

saveRDS(
  list(
    trial_counts_by_subject_condition = trial_counts_by_subject_condition,
    trial_level_descriptives = trial_level_descriptives,
    subject_level_means = subject_level_means,
    subject_level_descriptives = subject_level_descriptives,
    pairwise_contrasts = pairwise_contrasts,
    mixed_model_fixed_effects = mixed_model_stats,
    mixed_model_random_effects = random_effects_summary
  ),
  file = file.path(output_dir, "pilot_stats_for_power_analysis.rds")
)

message("Pilot statistics export complete.")
