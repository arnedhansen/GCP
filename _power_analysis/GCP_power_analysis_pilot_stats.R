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

## Collect singular-fit notices; printed and emitted as warnings at end of script.
pilot_warn_env <- new.env(parent = emptyenv())
pilot_warn_env$msgs <- character(0)
note_pilot_warning <- function(text) {
  pilot_warn_env$msgs <- c(pilot_warn_env$msgs, text)
}

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
PREFERRED_OUTPUT_DIR <- file.path(roots$data_root, "data", "pilot_stats")
FALLBACK_OUTPUT_DIR <- file.path(getwd(), "data", "pilot_stats")

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

resolve_merged_trials_csv <- function() {
  candidates <- c(
    file.path(roots$data_root, "data", "features", "GCP_merged_data_trials.csv"),
    file.path(getwd(), "data", "features", "GCP_merged_data_trials.csv")
  )
  hit <- candidates[file.exists(candidates)]
  if (length(hit) == 0) {
    return(NA_character_)
  }
  hit[1]
}
MERGED_TRIALS_FILE <- resolve_merged_trials_csv()

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
      sesoi_std = 0.05,
      external_lower_std = 0.05,
      external_point_std = 0.08
    ),
    PeakAmplitude = list(
      sesoi_std = -0.02,
      external_lower_std = -0.04,
      external_point_std = -0.06
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
      residual_sd_multiplier = c(1.50, 1.20, 1.35, 1.10),
      random_intercept_sd_multiplier = c(1.50, 1.20, 1.30, 1.10),
      random_slope_sd_multiplier = c(1.50, 1.20, 1.30, 1.10),
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

normalize_pilot_contrast_numbers <- function(cond) {
  v <- if (is.factor(cond)) as.character(cond) else cond
  num <- suppressWarnings(as.numeric(v))
  if (length(num) != length(cond)) {
    return(rep(NA_real_, length(cond)))
  }
  if (all(num %in% c(1L, 2L, 3L, 4L), na.rm = TRUE)) {
    return(c(25, 50, 75, 100)[as.integer(num)])
  }
  num
}

## Subject-level nonparametric bootstrap of residual SD (sigma) for the interaction lmer.
bootstrap_interaction_residual_sigma <- function(take, B = 500L, rng_seed = 123L) {
  set.seed(as.integer(rng_seed))
  subj_levels <- levels(take$Subject)
  n_sub <- length(subj_levels)
  sig <- rep(NA_real_, B)
  singular_boot <- 0L
  form <- gamma_power ~ contrast_num_c * microsaccade_c + (1 | Subject)
  for (b in seq_len(B)) {
    idx <- sample.int(n_sub, replace = TRUE)
    parts <- vector("list", length(idx))
    for (j in seq_along(idx)) {
      sb <- subj_levels[idx[j]]
      ch <- take[take$Subject == sb, , drop = FALSE]
      if (nrow(ch) == 0) {
        next
      }
      ch$Subject <- factor(rep(paste0(as.character(sb), ".__b", j), nrow(ch)))
      parts[[j]] <- ch
    }
    ok_parts <- parts[!vapply(parts, is.null, logical(1))]
    if (length(ok_parts) == 0) {
      next
    }
    take_b <- do.call(rbind, ok_parts)
    take_b$contrast_num_c <- as.numeric(scale(take_b$contrast_num, center = TRUE, scale = TRUE))
    take_b$microsaccade_c <- as.numeric(scale(take_b$microsaccade_raw, center = TRUE, scale = TRUE))
    if (nlevels(droplevels(take_b$Subject)) < 3L) {
      next
    }
    fit <- tryCatch(
      suppressMessages(lmer(form, data = take_b, REML = FALSE)),
      error = function(e) NULL
    )
    if (!is.null(fit)) {
      if (isTRUE(lme4::isSingular(fit, tol = 1e-4))) {
        singular_boot <- singular_boot + 1L
      }
      sg <- tryCatch(sigma(fit), error = function(e) NA_real_)
      if (is.finite(sg)) {
        sig[b] <- sg
      }
    }
  }
  ok <- is.finite(sig)
  rate <- mean(ok)
  n_ok <- sum(ok)
  if (n_ok < 10L) {
    stop(
      "Interaction residual bootstrap produced too few valid sigma draws (",
      n_ok, " < 10; success rate ", sprintf("%.2f", rate), ").",
      call. = FALSE
    )
  }
  qs <- stats::quantile(sig[ok], probs = c(0.025, 0.5, 0.975), names = FALSE)
  list(
    quantiles = qs,
    success_rate = rate,
    sigma_vector = sig,
    singular_replicates = singular_boot,
    B = B
  )
}

## Subject-level bootstrap of residual SD (sigma) for PeakFrequency / PeakAmplitude lmers (same
## formula as the final mm_freq / mm_power objects, including after singular refit).
bootstrap_subject_level_residual_sigma <- function(subject_level_means, formula_obj, B = 500L, rng_seed = 123L) {
  set.seed(as.integer(rng_seed))
  slm <- subject_level_means
  slm$Subject <- factor(slm$Subject)
  subj_levels <- levels(slm$Subject)
  n_sub <- length(subj_levels)
  sig <- rep(NA_real_, B)
  singular_boot <- 0L
  for (b in seq_len(B)) {
    idx <- sample.int(n_sub, replace = TRUE)
    parts <- vector("list", length(idx))
    for (j in seq_along(idx)) {
      sb <- subj_levels[idx[j]]
      ch <- slm[slm$Subject == sb, , drop = FALSE]
      if (nrow(ch) == 0) {
        next
      }
      ch$Subject <- factor(rep(paste0(as.character(sb), ".__b", j), nrow(ch)))
      parts[[j]] <- ch
    }
    ok_parts <- parts[!vapply(parts, is.null, logical(1))]
    if (length(ok_parts) == 0) {
      next
    }
    take_b <- do.call(rbind, ok_parts)
    take_b$contrast_num_c <- as.numeric(scale(take_b$contrast_num, center = TRUE, scale = TRUE))
    take_b$contrast_num_c2 <- take_b$contrast_num_c^2
    if (nlevels(droplevels(take_b$Subject)) < 3L) {
      next
    }
    fit <- tryCatch(
      suppressMessages(lmer(formula_obj, data = take_b, REML = FALSE)),
      error = function(e) NULL
    )
    if (!is.null(fit)) {
      if (isTRUE(lme4::isSingular(fit, tol = 1e-4))) {
        singular_boot <- singular_boot + 1L
      }
      sg <- tryCatch(sigma(fit), error = function(e) NA_real_)
      if (is.finite(sg)) {
        sig[b] <- sg
      }
    }
  }
  ok <- is.finite(sig)
  rate <- mean(ok)
  n_ok <- sum(ok)
  if (n_ok < 10L) {
    stop(
      "Subject-level residual bootstrap produced too few valid sigma draws (",
      n_ok, " < 10; success rate ", sprintf("%.2f", rate), ").",
      call. = FALSE
    )
  }
  qs <- stats::quantile(sig[ok], probs = c(0.025, 0.5, 0.975), names = FALSE)
  list(
    quantiles = qs,
    success_rate = rate,
    sigma_vector = sig,
    singular_replicates = singular_boot,
    B = B
  )
}

fit_and_export_interaction_pilot <- function(merged_path, output_dir) {
  if (is.na(merged_path) || length(merged_path) != 1L || !nzchar(merged_path) || !file.exists(merged_path)) {
    message("No GCP_merged_data_trials.csv found; skipping EEG–eye interaction pilot export.")
    return(NULL)
  }
  message("Loading merged trial table for interaction parameterization: ", merged_path)
  mt <- utils::read.csv(merged_path, stringsAsFactors = FALSE, check.names = TRUE)
  id_col <- if ("ID" %in% names(mt)) {
    "ID"
  } else if ("Subject" %in% names(mt)) {
    "Subject"
  } else {
    NA_character_
  }
  if (is.na(id_col) || !"Condition" %in% names(mt)) {
    message("Merged table missing ID/Subject or Condition; skipping interaction export.")
    return(NULL)
  }
  y_col <- intersect(c("GED_PeakAmplitude", "PeakAmplitude"), names(mt))[1]
  ms_col <- intersect(c("Gaze_PctMSRate", "Gaze_MSRate", "PctMSRate", "MSRate"), names(mt))[1]
  if (is.na(y_col) || is.na(ms_col)) {
    message("Merged table missing GED_PeakAmplitude or gaze microsaccade column; skipping interaction export.")
    return(NULL)
  }
  take <- mt[, c(id_col, "Condition", y_col, ms_col), drop = FALSE]
  names(take) <- c("Subject", "Condition", "gamma_power", "microsaccade_raw")
  take <- take[is.finite(take$gamma_power) & is.finite(take$microsaccade_raw), , drop = FALSE]
  take$Subject <- factor(take$Subject)
  take$contrast_num <- normalize_pilot_contrast_numbers(take$Condition)
  take <- take[is.finite(take$contrast_num), , drop = FALSE]
  if (nrow(take) < 50L || nlevels(take$Subject) < 3L) {
    message("Insufficient merged trials for interaction model (need >=50 trials, >=3 subjects). Skipping.")
    return(NULL)
  }
  take$contrast_num_c <- as.numeric(scale(take$contrast_num, center = TRUE, scale = TRUE))
  take$microsaccade_c <- as.numeric(scale(take$microsaccade_raw, center = TRUE, scale = TRUE))

  mm_ms <- tryCatch(
    suppressMessages(lmer(microsaccade_raw ~ contrast_num_c + (1 | Subject), data = take, REML = FALSE)),
    error = function(e) NULL
  )
  mm_int <- tryCatch(
    suppressMessages(lmer(gamma_power ~ contrast_num_c * microsaccade_c + (1 | Subject), data = take, REML = FALSE)),
    error = function(e) NULL
  )
  if (is.null(mm_int)) {
    message("Interaction mixed model failed to fit; skipping interaction export.")
    return(NULL)
  }
  if (isTRUE(lme4::isSingular(mm_int, tol = 1e-4))) {
    note_pilot_warning(
      "Interaction pilot model (gamma_power ~ contrast * microsaccade_c + (1|Subject), trial level): singular fit."
    )
  }
  nm_int <- "contrast_num_c:microsaccade_c"
  fix_i <- lme4::fixef(mm_int)
  if (!nm_int %in% names(fix_i)) {
    message("Interaction coefficient not in model summary; skipping interaction export.")
    return(NULL)
  }
  sm <- suppressMessages(summary(mm_int))
  se_ix <- tryCatch(
    as.numeric(sm$coefficients[nm_int, grep("Std\\.?\\s*Error", colnames(sm$coefficients), ignore.case = TRUE)[1]]),
    error = function(e) NA_real_
  )
  b_ix <- unname(fix_i[[nm_int]])
  b_con <- unname(fix_i[["contrast_num_c"]])
  b_ms <- unname(fix_i[["microsaccade_c"]])
  sigma_int <- sigma(mm_int)
  vc_int <- as.data.frame(VarCorr(mm_int))
  rid <- vc_int$grp == "Subject" & vc_int$var1 == "(Intercept)" & (is.na(vc_int$var2) | vc_int$var2 == "")
  ri_sd <- suppressWarnings(as.numeric(vc_int[rid, "sdcor"][1]))
  if (!is.finite(ri_sd)) {
    stop("Non-finite Subject random-intercept SD from interaction pilot model.", call. = FALSE)
  }

  if (is.null(mm_ms)) {
    stop("Microsaccade mixed model (microsaccade_raw ~ contrast_num_c + (1|Subject)) failed to fit.", call. = FALSE)
  }
  if (isTRUE(lme4::isSingular(mm_ms, tol = 1e-4))) {
    note_pilot_warning(
      "Microsaccade pilot model (microsaccade_raw ~ contrast_num_c + (1|Subject), trial level): singular fit."
    )
  }
  fix_ms <- lme4::fixef(mm_ms)
  b_ms_con <- unname(fix_ms[["contrast_num_c"]])
  sigma_ms <- sigma(mm_ms)
  vc_ms <- as.data.frame(VarCorr(mm_ms))
  rid_ms <- vc_ms$grp == "Subject" & vc_ms$var1 == "(Intercept)" & (is.na(vc_ms$var2) | vc_ms$var2 == "")
  ms_subj_sd <- suppressWarnings(as.numeric(vc_ms[rid_ms, "sdcor"][1]))
  if (!is.finite(ms_subj_sd)) {
    stop("Non-finite Subject random-intercept SD from microsaccade pilot model.", call. = FALSE)
  }

  ci_low_ix <- if (is.finite(se_ix)) b_ix - 1.96 * se_ix else NA_real_
  ci_high_ix <- if (is.finite(se_ix)) b_ix + 1.96 * se_ix else NA_real_

  INTERACTION_RESIDUAL_BOOTSTRAP_B <- 500L
  message(
    "Bootstrapping residual SD for interaction model (subject resampling, B=",
    INTERACTION_RESIDUAL_BOOTSTRAP_B, ") …"
  )
  bs_sig <- bootstrap_interaction_residual_sigma(
    take,
    B = INTERACTION_RESIDUAL_BOOTSTRAP_B,
    rng_seed = 123L
  )
  if (is.null(bs_sig$quantiles) || length(bs_sig$quantiles) != 3L) {
    stop(
      "Interaction residual bootstrap did not return three quantiles (success_rate=",
      sprintf("%.2f", bs_sig$success_rate), ").",
      call. = FALSE
    )
  }
  res_q025 <- as.numeric(bs_sig$quantiles[1])
  res_q50 <- as.numeric(bs_sig$quantiles[2])
  res_q975 <- as.numeric(bs_sig$quantiles[3])
  if (any(!is.finite(c(res_q025, res_q50, res_q975)))) {
    stop("Bootstrap residual quantiles are not all finite.", call. = FALSE)
  }
  message(
    "Bootstrap sigma quantiles (2.5%, 50%, 97.5%): ",
    sprintf("%.4f, %.4f, %.4f", res_q025, res_q50, res_q975),
    " | valid fits: ", sprintf("%.1f%%", 100 * bs_sig$success_rate)
  )
  n_sing_bs <- bs_sig$singular_replicates
  if (is.null(n_sing_bs)) {
    n_sing_bs <- 0L
  }
  if (n_sing_bs > 0L) {
    note_pilot_warning(sprintf(
      "Interaction residual bootstrap: %d of %d replicate lmer fits were singular (tol=1e-4).",
      as.integer(n_sing_bs),
      as.integer(INTERACTION_RESIDUAL_BOOTSTRAP_B)
    ))
  }

  ms_mean <- mean(take$microsaccade_raw, na.rm = TRUE)
  out <- data.frame(
    source_merged_csv = merged_path,
    microsaccade_column_used = ms_col,
    outcome_column_used = y_col,
    n_trials = nrow(take),
    n_subjects = nlevels(take$Subject),
    microsaccade_grand_mean = ms_mean,
    outcome_mean = mean(take$gamma_power, na.rm = TRUE),
    beta_contrast = b_con,
    beta_microsaccade_scaled = b_ms,
    beta_interaction = b_ix,
    beta_interaction_se = se_ix,
    beta_interaction_ci95_low = ci_low_ix,
    beta_interaction_ci95_high = ci_high_ix,
    random_intercept_sd = ri_sd,
    residual_sd_point = sigma_int,
    residual_sd_q025 = res_q025,
    residual_sd_q50 = res_q50,
    residual_sd_q975 = res_q975,
    residual_bootstrap_B = INTERACTION_RESIDUAL_BOOTSTRAP_B,
    residual_bootstrap_success_rate = bs_sig$success_rate,
    ms_beta_contrast_raw = b_ms_con,
    ms_random_intercept_sd = ms_subj_sd,
    ms_residual_sd = sigma_ms,
    stringsAsFactors = FALSE
  )
  utils::write.csv(
    out,
    file.path(output_dir, "pilot_interaction_power_parameters.csv"),
    row.names = FALSE
  )
  message("Wrote interaction pilot parameters: ", file.path(output_dir, "pilot_interaction_power_parameters.csv"))
  out
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
  note_pilot_warning(
    "PeakFrequency (subject-level means): initial model singular; refitting with (1 | Subject) only."
  )
  mm_freq <- lmer(
    PeakFrequency ~ contrast_num_c + contrast_num_c2 + (1 | Subject),
    data = subject_level_means,
    REML = FALSE
  )
  if (isSingular(mm_freq, tol = 1e-4)) {
    note_pilot_warning("PeakFrequency (subject-level means): refitted random-intercept-only model is still singular.")
  }
}

mm_power <- lmer(
  PeakAmplitude ~ contrast_num_c + contrast_num_c2 + (1 + contrast_num_c | Subject),
  data = subject_level_means,
  REML = FALSE
)
if (isSingular(mm_power, tol = 1e-4)) {
  note_pilot_warning(
    "PeakAmplitude (subject-level means): initial model singular; refitting with (1 | Subject) only."
  )
  mm_power <- lmer(
    PeakAmplitude ~ contrast_num_c + contrast_num_c2 + (1 | Subject),
    data = subject_level_means,
    REML = FALSE
  )
  if (isSingular(mm_power, tol = 1e-4)) {
    note_pilot_warning("PeakAmplitude (subject-level means): refitted random-intercept-only model is still singular.")
  }
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

SUBJECT_LEVEL_RESIDUAL_BOOTSTRAP_B <- 500L
message(
  "Bootstrapping residual SD for subject-level models (B=",
  SUBJECT_LEVEL_RESIDUAL_BOOTSTRAP_B, ", PeakFrequency) …"
)
bs_sigma_freq <- bootstrap_subject_level_residual_sigma(
  subject_level_means,
  stats::formula(mm_freq),
  B = SUBJECT_LEVEL_RESIDUAL_BOOTSTRAP_B,
  rng_seed = 123L
)
message(
  "  sigma q2.5 / q50 / q97.5: ",
  sprintf("%.4f, %.4f, %.4f", bs_sigma_freq$quantiles[1], bs_sigma_freq$quantiles[2], bs_sigma_freq$quantiles[3]),
  " | valid: ", sprintf("%.1f%%", 100 * bs_sigma_freq$success_rate)
)
if (bs_sigma_freq$singular_replicates > 0L) {
  note_pilot_warning(sprintf(
    "PeakFrequency subject-level residual bootstrap: %d of %d replicate lmer fits were singular (tol=1e-4).",
    as.integer(bs_sigma_freq$singular_replicates),
    as.integer(SUBJECT_LEVEL_RESIDUAL_BOOTSTRAP_B)
  ))
}
message(
  "Bootstrapping residual SD for subject-level models (B=",
  SUBJECT_LEVEL_RESIDUAL_BOOTSTRAP_B, ", PeakAmplitude) …"
)
bs_sigma_power <- bootstrap_subject_level_residual_sigma(
  subject_level_means,
  stats::formula(mm_power),
  B = SUBJECT_LEVEL_RESIDUAL_BOOTSTRAP_B,
  rng_seed = 124L
)
message(
  "  sigma q2.5 / q50 / q97.5: ",
  sprintf("%.4f, %.4f, %.4f", bs_sigma_power$quantiles[1], bs_sigma_power$quantiles[2], bs_sigma_power$quantiles[3]),
  " | valid: ", sprintf("%.1f%%", 100 * bs_sigma_power$success_rate)
)
if (bs_sigma_power$singular_replicates > 0L) {
  note_pilot_warning(sprintf(
    "PeakAmplitude subject-level residual bootstrap: %d of %d replicate lmer fits were singular (tol=1e-4).",
    as.integer(bs_sigma_power$singular_replicates),
    as.integer(SUBJECT_LEVEL_RESIDUAL_BOOTSTRAP_B)
  ))
}
pilot_subject_level_residual_bootstrap <- data.frame(
  outcome = c("PeakFrequency", "PeakAmplitude"),
  residual_sd_q025 = c(bs_sigma_freq$quantiles[1], bs_sigma_power$quantiles[1]),
  residual_sd_q50 = c(bs_sigma_freq$quantiles[2], bs_sigma_power$quantiles[2]),
  residual_sd_q975 = c(bs_sigma_freq$quantiles[3], bs_sigma_power$quantiles[3]),
  residual_bootstrap_B = SUBJECT_LEVEL_RESIDUAL_BOOTSTRAP_B,
  residual_bootstrap_success_rate = c(bs_sigma_freq$success_rate, bs_sigma_power$success_rate),
  singular_replicates = c(bs_sigma_freq$singular_replicates, bs_sigma_power$singular_replicates),
  lmer_formula = c(
    paste(trimws(deparse(stats::formula(mm_freq))), collapse = " "),
    paste(trimws(deparse(stats::formula(mm_power))), collapse = " ")
  ),
  stringsAsFactors = FALSE
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
write.csv(
  triangulated_manifest,
  file.path(output_dir, "power_parameter_manifest.csv"),
  row.names = FALSE
)
write.csv(
  pilot_subject_level_residual_bootstrap,
  file.path(output_dir, "pilot_subject_level_residual_sigma_bootstrap.csv"),
  row.names = FALSE
)
message(
  "Wrote: ",
  file.path(output_dir, "pilot_subject_level_residual_sigma_bootstrap.csv")
)

interaction_pilot_params <- fit_and_export_interaction_pilot(MERGED_TRIALS_FILE, output_dir)

saveRDS(
  list(
    trial_counts_by_subject_condition = trial_counts_by_subject_condition,
    trial_level_descriptives = trial_level_descriptives,
    subject_level_means = subject_level_means,
    subject_level_descriptives = subject_level_descriptives,
    pairwise_contrasts = pairwise_contrasts,
    mixed_model_fixed_effects = mixed_model_stats,
    mixed_model_random_effects = random_effects_summary,
    triangulated_manifest = triangulated_manifest,
    pilot_subject_level_residual_bootstrap = pilot_subject_level_residual_bootstrap,
    interaction_pilot_params = interaction_pilot_params
  ),
  file = file.path(output_dir, "pilot_stats_for_power_analysis.rds")
)

message("Pilot statistics export complete.")

if (length(pilot_warn_env$msgs) > 0) {
  message("\n========== PILOT WARNINGS: singular or simplified mixed-model fits ==========")
  for (t in pilot_warn_env$msgs) {
    message("  • ", t)
  }
  message("============================================================================")
  for (t in pilot_warn_env$msgs) {
    base::warning(t, call. = FALSE)
  }
} else {
  message("No singular pilot mixed-model fits were flagged (isSingular, tol=1e-4).")
}
