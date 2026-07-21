## Shared helpers for GCP simulation-based power analysis.
## Each simulated/fitted dataset has exactly one row per subject per contrast condition.

SUBJECT_CONDITIONS_PER_SUBJECT <- 4L

resolve_gcp_root <- function() {
  candidates <- if (.Platform$OS.type == "windows") {
    c("W:/Students/Arne/GCP", getwd())
  } else {
    c(
      "/Volumes/g_psyplafor_methlab$/Students/Arne/GCP",
      file.path(getwd()),
      file.path(dirname(getwd()))
    )
  }
  probe <- file.path(candidates, "data", "pilot_stats", "pilot_subject_level_variance_partition_summary.csv")
  hit <- candidates[file.exists(probe)]
  if (length(hit) > 0) {
    return(hit[1])
  }
  candidates[1]
}

load_residual_sd_scenarios <- function(gcp_root, pilot_outcome, fallback_levels, model_scope = NULL) {
  partition <- tryCatch(
    load_variance_partition_summary(gcp_root, pilot_outcome, model_scope),
    error = function(e) NULL
  )
  if (!is.null(partition)) {
    levels <- partition$component_quantiles$residual_sd
    if (all(is.finite(levels)) && all(levels > 0)) {
      baseline <- levels[2L]
      return(list(
        residual_sd_levels = levels,
        baseline_residual_sd = baseline,
        residual_multipliers = levels / baseline,
        source = partition$source_path,
        data_source = partition$data_source
      ))
    }
  }
  baseline <- fallback_levels[2L]
  list(
    residual_sd_levels = fallback_levels,
    baseline_residual_sd = baseline,
    residual_multipliers = fallback_levels / baseline,
    source = "hardcoded_fallback"
  )
}

safe_random_slope_sd <- function(model_obj, term_name) {
  vc <- as.data.frame(lme4::VarCorr(model_obj))
  hit <- vc$grp == "Subject" & vc$var1 == term_name & (is.na(vc$var2) | vc$var2 == "")
  if (sum(hit) == 0) {
    return(0)
  }
  as.numeric(vc$sdcor[hit][1])
}

load_pilot_random_effects <- function(gcp_root, pilot_outcome, defaults, model_scope = NULL) {
  re_path <- file.path(
    gcp_root, "data", "pilot_stats",
    "pilot_mixed_model_random_effects_gamma_freq_power.csv"
  )
  fe_path <- file.path(
    gcp_root, "data", "pilot_stats",
    "pilot_mixed_model_fixed_effects_gamma_freq_power.csv"
  )
  slm_path <- file.path(
    gcp_root, "data", "pilot_stats",
    "pilot_subject_level_means_gamma_freq_power.csv"
  )

  out <- defaults
  out$source <- "hardcoded_fallback"

  if (!file.exists(re_path) || !file.exists(fe_path)) {
    return(out)
  }

  re <- utils::read.csv(re_path, stringsAsFactors = FALSE)
  fe <- utils::read.csv(fe_path, stringsAsFactors = FALSE)
  if (!is.null(model_scope) && "model_scope" %in% names(re)) {
    re <- re[re$model_scope == model_scope, , drop = FALSE]
    fe <- fe[fe$model_scope == model_scope, , drop = FALSE]
  }
  re <- re[re$outcome == pilot_outcome, , drop = FALSE]
  fe <- fe[fe$outcome == pilot_outcome, , drop = FALSE]

  ri_hit <- re$grp == "Subject" & re$var1 == "(Intercept)" & (is.na(re$var2) | re$var2 == "")
  if (any(ri_hit)) {
    out$random_intercept_sd <- as.numeric(re$sdcor[ri_hit][1])
  }

  rs_hit <- re$grp == "Subject" & re$var1 == "contrast_num_c" & (is.na(re$var2) | re$var2 == "")
  if (any(rs_hit)) {
    out$random_slope_sd <- as.numeric(re$sdcor[rs_hit][1])
  }

  rqs_hit <- re$grp == "Subject" & re$var1 == "contrast_num_c2" & (is.na(re$var2) | re$var2 == "")
  if (any(rqs_hit) && !is.null(defaults$random_quadratic_slope_sd)) {
    out$random_quadratic_slope_sd <- as.numeric(re$sdcor[rqs_hit][1])
  }

  intercept_fe <- fe$term == "(Intercept)"
  if (any(intercept_fe)) {
    out$outcome_mean <- as.numeric(fe$estimate[intercept_fe][1])
  }

  lin_fe <- fe$term == "contrast_num_c"
  if (any(lin_fe) && !is.null(defaults$linear_nuisance_beta)) {
    out$linear_nuisance_beta <- as.numeric(fe$estimate[lin_fe][1])
  }

  if (file.exists(slm_path)) {
    slm <- utils::read.csv(slm_path, stringsAsFactors = FALSE)
    value_col <- if (pilot_outcome == "PeakFrequency") "PeakFrequency" else "PeakAmplitude"
    if (value_col %in% names(slm)) {
      out$outcome_mean <- mean(slm[[value_col]], na.rm = TRUE)
    }
  }

  out$source <- "pilot_stats_csv"
  out
}

extract_variance_partition_from_fit <- function(fit) {
  vc <- as.data.frame(lme4::VarCorr(fit))
  ri_hit <- vc$grp == "Subject" & vc$var1 == "(Intercept)" &
    (is.na(vc$var2) | vc$var2 == "")
  random_intercept_sd <- if (any(ri_hit)) as.numeric(vc$sdcor[ri_hit][1]) else 0
  list(
    residual_sd = stats::sigma(fit),
    random_intercept_sd = random_intercept_sd,
    random_slope_sd = safe_random_slope_sd(fit, "contrast_num_c"),
    random_quadratic_slope_sd = safe_random_slope_sd(fit, "contrast_num_c2")
  )
}

infer_random_effect_structure <- function(formula_label) {
  form <- as.character(formula_label)
  includes_random_intercept <- grepl("\\(1", form)
  includes_random_slope <- grepl("\\(\\s*1\\s*\\+\\s*contrast_num_c", form) ||
    grepl("0\\s*\\+\\s*contrast_num_c\\s*\\|", form)
  includes_random_quadratic_slope <- grepl("contrast_num_c2", form) && (
    grepl("\\+\\s*contrast_num_c2", form) || grepl("0\\s*\\+\\s*contrast_num_c2", form)
  )
  list(
    includes_random_intercept = includes_random_intercept,
    includes_random_slope = includes_random_slope,
    includes_random_quadratic_slope = includes_random_quadratic_slope
  )
}

component_quantile <- function(x, probs = c(0.025, 0.5, 0.975)) {
  x <- x[is.finite(x)]
  if (length(x) < 10L) {
    return(rep(NA_real_, length(probs)))
  }
  as.numeric(stats::quantile(x, probs = probs, names = FALSE))
}

`%||%` <- function(x, y) if (is.null(x)) y else x

bootstrap_subject_level_variance_partition <- function(
    subject_level_means,
    formula_obj,
    B = 500L,
    rng_seed = 123L) {
  set.seed(as.integer(rng_seed))
  slm <- subject_level_means
  slm$Subject <- factor(slm$Subject)
  subj_levels <- levels(slm$Subject)
  n_sub <- length(subj_levels)
  partitions <- vector("list", B)
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
      suppressMessages(lme4::lmer(
        formula_obj,
        data = take_b,
        REML = FALSE,
        control = lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
      )),
      error = function(e) NULL
    )
    if (!is.null(fit)) {
      if (isTRUE(lme4::isSingular(fit, tol = 1e-4))) {
        singular_boot <- singular_boot + 1L
      }
      partitions[[b]] <- extract_variance_partition_from_fit(fit)
    }
  }
  ok_idx <- vapply(partitions, Negate(is.null), logical(1))
  if (sum(ok_idx) < 10L) {
    stop(
      "Variance-partition bootstrap produced too few valid fits (",
      sum(ok_idx), " < 10)."
    )
  }
  mat <- do.call(rbind, lapply(partitions[ok_idx], function(x) {
    data.frame(
      residual_sd = x$residual_sd,
      random_intercept_sd = x$random_intercept_sd,
      random_slope_sd = x$random_slope_sd,
      random_quadratic_slope_sd = x$random_quadratic_slope_sd,
      stringsAsFactors = FALSE
    )
  }))
  list(
    draws = as.data.frame(mat),
    success_rate = mean(ok_idx),
    singular_replicates = singular_boot,
    B = B
  )
}

summarize_variance_partition_bootstrap <- function(
    bs_obj,
    outcome,
    model_scope,
    data_source,
    pilot_n_subjects,
    lmer_formula,
    outcome_mean = NA_real_,
    linear_nuisance_beta = NA_real_) {
  cmp_names <- c(
    "residual_sd", "random_intercept_sd", "random_slope_sd", "random_quadratic_slope_sd"
  )
  meta_cols <- c(
    "data_source", "pilot_n_subjects", "lmer_formula", "outcome_mean",
    "linear_nuisance_beta", "bootstrap_B", "bootstrap_success_rate", "singular_replicates"
  )
  cmp_rows <- lapply(cmp_names, function(cmp) {
    x <- bs_obj$draws[[cmp]]
    q <- component_quantile(x)
    row <- data.frame(
      record_type = "component",
      outcome = outcome,
      model_scope = model_scope,
      component = cmp,
      q025 = q[1], q50 = q[2], q975 = q[3],
      stringsAsFactors = FALSE
    )
    for (nm in meta_cols) {
      row[[nm]] <- NA
    }
    row
  })
  meta_row <- data.frame(
    record_type = "meta",
    outcome = outcome,
    model_scope = model_scope,
    component = NA_character_,
    q025 = NA_real_, q50 = NA_real_, q975 = NA_real_,
    data_source = data_source,
    pilot_n_subjects = pilot_n_subjects,
    lmer_formula = lmer_formula,
    outcome_mean = outcome_mean,
    linear_nuisance_beta = linear_nuisance_beta,
    bootstrap_B = bs_obj$B,
    bootstrap_success_rate = bs_obj$success_rate,
    singular_replicates = bs_obj$singular_replicates,
    stringsAsFactors = FALSE
  )
  do.call(rbind, c(list(meta_row), cmp_rows))
}

bootstrap_subject_level_residual_sigma <- function(
    subject_level_means, formula_obj, B = 500L, rng_seed = 123L) {
  bs <- bootstrap_subject_level_variance_partition(
    subject_level_means, formula_obj, B = B, rng_seed = rng_seed
  )
  qs <- component_quantile(bs$draws$residual_sd)
  list(
    quantiles = qs,
    success_rate = bs$success_rate,
    sigma_vector = bs$draws$residual_sd,
    singular_replicates = bs$singular_replicates,
    B = B
  )
}

load_variance_partition_summary <- function(gcp_root, pilot_outcome, model_scope) {
  summary_path <- file.path(
    gcp_root, "data", "pilot_stats",
    "pilot_subject_level_variance_partition_summary.csv"
  )
  if (!file.exists(summary_path)) {
    stop("Variance partition summary not found: ", summary_path)
  }
  sm <- utils::read.csv(summary_path, stringsAsFactors = FALSE)
  meta <- sm[sm$outcome == pilot_outcome & sm$model_scope == model_scope & sm$record_type == "meta", , drop = FALSE]
  cmp <- sm[sm$outcome == pilot_outcome & sm$model_scope == model_scope & sm$record_type == "component", , drop = FALSE]
  if (nrow(meta) != 1L) {
    stop("Expected one meta row for ", pilot_outcome, " / ", model_scope)
  }

  qmat <- function(component) {
    row <- cmp[cmp$component == component, , drop = FALSE]
    if (nrow(row) != 1L) {
      return(c(NA_real_, NA_real_, NA_real_))
    }
    c(row$q025, row$q50, row$q975)
  }

  re_struct <- infer_random_effect_structure(meta$lmer_formula[1])
  list(
    source_path = summary_path,
    data_source = meta$data_source[1],
    pilot_n_subjects = as.integer(meta$pilot_n_subjects[1]),
    lmer_formula = meta$lmer_formula[1],
    bootstrap_B = as.integer(meta$bootstrap_B[1]),
    singular_replicates = as.integer(meta$singular_replicates[1]),
    includes_random_intercept = re_struct$includes_random_intercept,
    includes_random_slope = re_struct$includes_random_slope,
    includes_random_quadratic_slope = re_struct$includes_random_quadratic_slope,
    component_quantiles = list(
      residual_sd = qmat("residual_sd"),
      random_intercept_sd = qmat("random_intercept_sd"),
      random_slope_sd = qmat("random_slope_sd"),
      random_quadratic_slope_sd = qmat("random_quadratic_slope_sd")
    ),
    point_estimates = list(
      outcome_mean = as.numeric(meta$outcome_mean[1]),
      linear_nuisance_beta = as.numeric(meta$linear_nuisance_beta[1])
    )
  )
}

pick_quantile <- function(qvec, level = c("low", "median", "high")) {
  level <- match.arg(level)
  idx <- switch(level, low = 1L, median = 2L, high = 3L)
  as.numeric(qvec[idx])
}

load_manifest_beta <- function(
    gcp_root,
    pilot_outcome,
    target_term,
    scenario_roles = c("sesoi", "pilot_secondary", "base")) {
  manifest_path <- file.path(gcp_root, "data", "pilot_stats", "power_parameter_manifest.csv")
  if (!file.exists(manifest_path)) {
    return(numeric(0))
  }
  manifest <- utils::read.csv(manifest_path, stringsAsFactors = FALSE)
  vals <- numeric(0)
  for (role in scenario_roles) {
    hit <- manifest$outcome == pilot_outcome &
      manifest$target_term == target_term &
      manifest$scenario_role == role
    if (sum(hit) == 1L && is.finite(manifest$beta_raw[hit][1])) {
      vals <- c(vals, as.numeric(manifest$beta_raw[hit][1]))
    }
  }
  vals
}

conservative_simulation_beta <- function(
    pilot_beta,
    manifest_betas,
    direction = c("positive", "negative")) {
  direction <- match.arg(direction)
  candidates <- c(pilot_beta, manifest_betas)
  candidates <- candidates[is.finite(candidates)]
  if (length(candidates) == 0L) {
    stop("No finite candidate betas for conservative simulation effect.")
  }
  if (direction == "positive") {
    pos <- candidates[candidates > 0]
    if (length(pos) > 0L) {
      return(min(pos))
    }
  }
  if (direction == "negative") {
    neg <- candidates[candidates < 0]
    if (length(neg) > 0L) {
      return(max(neg))
    }
  }
  candidates[which.min(abs(candidates))]
}

load_conservative_simulation_beta <- function(
    gcp_root,
    pilot_outcome,
    target_term,
    direction = c("positive", "negative"),
    pilot_fallback = NA_real_) {
  pilot_beta <- load_pilot_fixed_effect(gcp_root, pilot_outcome, target_term, pilot_fallback)
  manifest_betas <- load_manifest_beta(gcp_root, pilot_outcome, target_term)
  conservative_simulation_beta(pilot_beta, manifest_betas, direction)
}

cohen_d_to_contrast_beta <- function(cohen_d, outcome_sd) {
  as.numeric(cohen_d) * as.numeric(outcome_sd)
}

load_outcome_sd <- function(gcp_root, pilot_outcome) {
  manifest_path <- file.path(gcp_root, "data", "pilot_stats", "power_parameter_manifest.csv")
  if (file.exists(manifest_path)) {
    manifest <- utils::read.csv(manifest_path, stringsAsFactors = FALSE)
    hit <- manifest$outcome == pilot_outcome & manifest$scenario_role == "base"
    if (sum(hit) >= 1L && is.finite(manifest$outcome_sd[hit][1])) {
      return(as.numeric(manifest$outcome_sd[hit][1]))
    }
  }
  slm_path <- file.path(gcp_root, "data", "pilot_stats", "pilot_subject_level_means_gamma_freq_power.csv")
  if (file.exists(slm_path)) {
    slm <- utils::read.csv(slm_path, stringsAsFactors = FALSE)
    col <- if (pilot_outcome == "PeakFrequency") "PeakFrequency" else "PeakAmplitude"
    if (col %in% names(slm)) {
      return(stats::sd(slm[[col]], na.rm = TRUE))
    }
  }
  if (pilot_outcome == "PeakFrequency") 6.656 else 1.925
}

load_pilot_fixed_effect <- function(gcp_root, pilot_outcome, term, fallback = NA_real_) {
  fe_path <- file.path(gcp_root, "data", "pilot_stats", "pilot_mixed_model_fixed_effects_gamma_freq_power.csv")
  if (!file.exists(fe_path)) {
    return(fallback)
  }
  fe <- utils::read.csv(fe_path, stringsAsFactors = FALSE)
  scope <- if (term == "contrast_num_c") "linear" else "quadratic"
  hit <- fe$outcome == pilot_outcome & fe$model_scope == scope & fe$term == term
  if (sum(hit) == 1L && is.finite(fe$estimate[hit][1])) {
    return(as.numeric(fe$estimate[hit][1]))
  }
  fallback
}

build_pessimistic_scenario <- function(partition) {
  q <- partition$component_quantiles
  use_rs <- partition$includes_random_slope
  use_rqs <- partition$includes_random_quadratic_slope
  data.frame(
    scenario_label = "pessimistic",
    scenario_display = "Pessimistic",
    varied_component = "joint_high_variance",
    random_intercept_sd = pick_quantile(q$random_intercept_sd, "high"),
    random_slope_sd = if (use_rs) pick_quantile(q$random_slope_sd, "high") else 0,
    random_quadratic_slope_sd = if (use_rqs) pick_quantile(q$random_quadratic_slope_sd, "high") else 0,
    residual_sd = pick_quantile(q$residual_sd, "high"),
    stringsAsFactors = FALSE
  )
}

build_power_scenarios <- function(partition, sensitivity_axis = c("residual")) {
  sensitivity_axis <- match.arg(sensitivity_axis)
  q <- partition$component_quantiles
  use_rs <- partition$includes_random_slope
  use_rqs <- partition$includes_random_quadratic_slope

  ri_med <- pick_quantile(q$random_intercept_sd, "median")
  rs_med <- if (use_rs) pick_quantile(q$random_slope_sd, "median") else 0
  rqs_med <- if (use_rqs) pick_quantile(q$random_quadratic_slope_sd, "median") else 0
  res_med <- pick_quantile(q$residual_sd, "median")

  if (sensitivity_axis == "residual") {
    levels <- c("low", "median", "high")
    res_vals <- c(
      pick_quantile(q$residual_sd, "low"),
      pick_quantile(q$residual_sd, "median"),
      pick_quantile(q$residual_sd, "high")
    )
    scenario_df <- data.frame(
      scenario_label = paste0("res_", sprintf("%.2f", res_vals / res_med)),
      scenario_display = sprintf("Residual %.3f", res_vals),
      varied_component = "residual_only",
      random_intercept_sd = ri_med,
      random_slope_sd = rs_med,
      random_quadratic_slope_sd = rqs_med,
      residual_sd = res_vals,
      stringsAsFactors = FALSE
    )
    return(scenario_df)
  }

  stop("Only residual sensitivity axis is supported.")
}

format_re_scenario_display <- function(ri, rs, rqs, use_rs, use_rqs) {
  parts <- sprintf("RI %.2f", ri)
  if (use_rs) {
    parts <- c(parts, sprintf("RS %.3f", rs))
  }
  if (use_rqs) {
    parts <- c(parts, sprintf("RQS %.3f", rqs))
  }
  paste(parts, collapse = " | ")
}

fit_lmer_with_fallbacks <- function(formula_list, data) {
  assert_one_row_per_subject_condition(data)
  for (i in seq_along(formula_list)) {
    form <- formula_list[[i]]
    fit <- tryCatch(
      suppressMessages(lme4::lmer(
        form,
        data = data,
        REML = FALSE,
        control = lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
      )),
      error = function(e) NULL
    )
    if (!is.null(fit)) {
      singular <- isTRUE(lme4::isSingular(fit, tol = 1e-4))
      if (!singular || i == length(formula_list)) {
        return(list(
          fit = fit,
          formula_label = paste(trimws(deparse(form)), collapse = " "),
          singular = singular
        ))
      }
    }
  }
  NULL
}

make_subject_level_design <- function(n_subjects, contrast_levels) {
  # One row per subject per contrast condition (no repeated trials).
  Subject <- factor(rep(seq_len(n_subjects), each = length(contrast_levels)))
  contrast <- factor(
    rep(contrast_levels, times = n_subjects),
    levels = contrast_levels,
    ordered = TRUE
  )
  contrast_num <- as.numeric(as.character(contrast))
  contrast_num_c <- as.numeric(scale(contrast_num, center = TRUE, scale = TRUE))
  out <- data.frame(
    Subject = Subject,
    Contrast = contrast,
    contrast_num = contrast_num,
    contrast_num_c = contrast_num_c
  )
  out$contrast_num_c2 <- out$contrast_num_c^2
  out
}

assert_one_row_per_subject_condition <- function(
    dat,
    subject_col = "Subject",
    contrast_col = "contrast_num") {
  if (!subject_col %in% names(dat) || !contrast_col %in% names(dat)) {
    stop("Data must contain Subject and contrast identifiers.")
  }
  counts <- table(dat[[subject_col]], dat[[contrast_col]])
  if (any(counts > 1L)) {
    stop(
      "Expected exactly one observation per subject per condition; ",
      "found up to ", max(counts), " rows in a subject-condition cell."
    )
  }
  invisible(TRUE)
}

minimum_subject_level_rows <- function(n_subjects, conditions_per_subject = SUBJECT_CONDITIONS_PER_SUBJECT) {
  as.integer(n_subjects) * as.integer(conditions_per_subject)
}

subject_level_fit_eligible <- function(dat, min_subjects = 3L) {
  dat <- droplevels(dat)
  n_sub <- nlevels(dat$Subject)
  if (n_sub < min_subjects) {
    return(FALSE)
  }
  if (nrow(dat) != minimum_subject_level_rows(n_sub)) {
    return(FALSE)
  }
  tryCatch({
    assert_one_row_per_subject_condition(dat)
    TRUE
  }, error = function(e) FALSE)
}

build_subject_level_formula_list <- function(outcome_col, model_scope = c("linear", "quadratic")) {
  model_scope <- match.arg(model_scope)
  if (model_scope == "linear") {
    return(list(
      stats::as.formula(paste(outcome_col, "~ contrast_num_c + (1 + contrast_num_c || Subject)")),
      stats::as.formula(paste(outcome_col, "~ contrast_num_c + (1 + contrast_num_c | Subject)")),
      stats::as.formula(paste(outcome_col, "~ contrast_num_c + (1 | Subject)"))
    ))
  }
  list(
    stats::as.formula(paste(
      outcome_col,
      "~ contrast_num_c + contrast_num_c2 + (1 + contrast_num_c + contrast_num_c2 || Subject)"
    )),
    stats::as.formula(paste(
      outcome_col,
      "~ contrast_num_c + contrast_num_c2 + (1 + contrast_num_c || Subject)"
    )),
    stats::as.formula(paste(outcome_col, "~ contrast_num_c + contrast_num_c2 + (1 | Subject)"))
  )
}

coef_table_from_fit <- function(fit) {
  if (requireNamespace("lmerTest", quietly = TRUE)) {
    sm <- tryCatch(lmerTest::summary(fit), error = function(e) NULL)
    if (!is.null(sm) && !is.null(sm$coefficients) && nrow(sm$coefficients) > 0) {
      return(sm$coefficients)
    }
  }
  tryCatch(summary(fit)$coefficients, error = function(e) NULL)
}

detection_from_fit <- function(fit, term, detection_alpha = 0.05) {
  coef_tbl <- coef_table_from_fit(fit)
  if (is.null(coef_tbl) || !term %in% rownames(coef_tbl)) {
    return(list(detection_success = FALSE, est = NA_real_, valid = FALSE))
  }
  est <- as.numeric(coef_tbl[term, "Estimate"])
  est_se <- as.numeric(coef_tbl[term, "Std. Error"])
  p_col <- grep("^Pr\\(>\\|", colnames(coef_tbl), value = TRUE)
  p_val <- if (length(p_col) > 0) {
    as.numeric(coef_tbl[term, p_col[1]])
  } else {
    stats::pnorm(-abs(est / est_se)) * 2
  }
  list(
    detection_success = is.finite(p_val) && p_val < detection_alpha,
    est = est,
    valid = TRUE
  )
}

power_binomial_se <- function(successes, n) {
  if (n <= 0) {
    return(NA_real_)
  }
  p <- successes / n
  sqrt(p * (1 - p) / n)
}

n_for_target_power <- function(power_df, target = 0.90) {
  df <- power_df[order(power_df$n_subjects), , drop = FALSE]
  hit <- df[df$power >= target, "n_subjects", drop = TRUE]
  if (length(hit) > 0) {
    return(as.integer(min(hit)))
  }
  NA_integer_
}

save_power_figures <- function(
    power_df,
    scenario_levels,
    subject_breaks,
    strict_power_target,
    figure_output_dir,
    output_prefix,
    plot_title,
    y_legend_title = "Residuals SD") {
  heat_low_color <- "#8E0F1F"
  heat_high_color <- "#0B6E3A"
  plot_df <- power_df
  plot_df$scenario_display <- factor(plot_df$scenario_display, levels = scenario_levels)
  n_unique_n <- length(unique(plot_df$n_subjects))

  curve_plot <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = .data$n_subjects,
      y = .data$power,
      color = .data$scenario_display,
      group = .data$scenario_display
    )
  )
  if (n_unique_n > 1L) {
    curve_plot <- curve_plot + ggplot2::geom_line(linewidth = 0.9, linetype = "dotted")
  }
  curve_plot <- curve_plot +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
      linewidth = 0.6,
      width = if (n_unique_n > 1L) 1.6 else 0,
      alpha = 0.70
    ) +
    ggplot2::geom_point(size = 2.8) +
    ggplot2::geom_hline(
      yintercept = strict_power_target,
      linetype = "dashed", color = "grey45", linewidth = 0.7
    ) +
    ggplot2::scale_color_manual(
      values = stats::setNames(c("#1B9E77", "#D95F02", "#7570B3"), scenario_levels),
      breaks = scenario_levels,
      labels = scenario_levels
    ) +
    ggplot2::scale_x_continuous(breaks = sort(unique(subject_breaks))) +
    ggplot2::scale_y_continuous(
      labels = scales::percent_format(accuracy = 1),
      limits = c(0, 1),
      breaks = c(0, 0.25, 0.50, 0.75, 0.90, 1)
    ) +
    ggplot2::labs(x = "Subjects", y = "Power", color = y_legend_title) +
    ggplot2::theme_classic(base_size = 18, base_family = "Arial") +
    ggplot2::theme(
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.grid = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(color = "black", linewidth = 1.2),
      axis.ticks = ggplot2::element_line(color = "black", linewidth = 1.0),
      axis.ticks.length = grid::unit(0.25, "cm"),
      axis.title = ggplot2::element_text(size = 20),
      axis.text = ggplot2::element_text(size = 16, color = "black"),
      legend.position = "right",
      legend.title = ggplot2::element_text(size = 15),
      legend.text = ggplot2::element_text(size = 11)
    )

  grDevices::png(
    file = file.path(figure_output_dir, paste0(output_prefix, ".png")),
    width = 2200, height = 1400, res = 300
  )
  print(curve_plot)
  grDevices::dev.off()

  heatmap_plot <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = factor(.data$n_subjects),
      y = .data$scenario_display,
      fill = .data$power
    )
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 1.1) +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.2f", .data$power)),
      color = "white", size = 3.52, fontface = "bold", family = "Arial"
    ) +
    ggplot2::scale_fill_gradientn(
      colours = c(heat_low_color, "#F46D43", "#FEE08B", "#66BD63", heat_high_color),
      values = c(0.00, 0.60, 0.79, 0.80, 1.00),
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.2),
      labels = sprintf("%.2f", seq(0, 1, by = 0.2))
    ) +
    ggplot2::coord_fixed() +
    ggplot2::labs(
      x = "Subjects", y = y_legend_title, fill = "Power"
    ) +
    ggplot2::theme_minimal(base_size = 16, base_family = "Arial") +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      axis.title = ggplot2::element_text(size = 13),
      axis.text = ggplot2::element_text(size = 10),
      legend.position = "right",
      legend.title = ggplot2::element_text(size = 11),
      legend.text = ggplot2::element_text(size = 10),
      panel.spacing = grid::unit(0.9, "lines")
    )

  grDevices::png(
    file = file.path(figure_output_dir, paste0(output_prefix, "_heatmap.png")),
    width = 2500, height = 1400, res = 300
  )
  print(heatmap_plot)
  grDevices::dev.off()
}
