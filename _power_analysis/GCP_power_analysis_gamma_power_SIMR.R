## GCP Simulation for Power Analysis on GAMMA POWER using SIMR
library(lme4)
library(simr)
library(ggplot2)

set.seed(123)

## Exploration vs reporting precision.
NSIM_EXPLORATORY <- 1000
NSIM_REPORTING <- 5000
NSIM_ACTIVE <- NSIM_EXPLORATORY
ALPHA <- 0.05

## Planned RR design assumptions.
CONTRAST_LEVELS <- c("25", "50", "75", "100")
N_SUBJECT_BASE <- 50
TRIALS_PER_CONDITION <- 160
SUBJECT_BREAKS <- seq(10, 100, by = 10)

## Effect assumptions for power curvature:
## - Uses plausible negative curvature values without over-anchoring to pilot.
effect_scenarios <- data.frame(
  effect_label = c("SESOI", "Conservative", "Optimistic"),
  beta_quadratic = c(-0.10, -0.20, -0.30)
)
random_slope_sd_scenarios <- c(0.05, 0.10, 0.15)

make_design <- function(n_subjects, trials_per_condition = TRIALS_PER_CONDITION) {
  Subject <- factor(rep(seq_len(n_subjects), each = length(CONTRAST_LEVELS) * trials_per_condition))
  contrast <- factor(
    rep(rep(CONTRAST_LEVELS, each = trials_per_condition), times = n_subjects),
    levels = CONTRAST_LEVELS,
    ordered = TRUE
  )
  contrast_num <- as.numeric(as.character(contrast))
  contrast_num_c <- as.numeric(scale(contrast_num, center = TRUE, scale = TRUE))
  contrast_num_c2 <- contrast_num_c^2
  data.frame(
    Subject = Subject,
    contrast = contrast,
    contrast_num = contrast_num,
    contrast_num_c = contrast_num_c,
    contrast_num_c2 = contrast_num_c2
  )
}

extract_warning_diagnostics <- function(warnings_vec) {
  if (length(warnings_vec) == 0) {
    return(list(n_warnings = 0L, convergence_rate = 0, singular_rate = 0))
  }
  convergence_hits <- grepl("converg|failed to converge|unable to evaluate", warnings_vec, ignore.case = TRUE)
  singular_hits <- grepl("singular", warnings_vec, ignore.case = TRUE)
  list(
    n_warnings = length(warnings_vec),
    convergence_rate = mean(convergence_hits),
    singular_rate = mean(singular_hits)
  )
}

run_power_with_diagnostics <- function(model_obj, test_term, nsim, along, breaks) {
  warn_bucket <- character(0)
  msg_bucket <- character(0)
  with_warning_capture <- function(expr) {
    withCallingHandlers(
      expr,
      warning = function(w) {
        warn_bucket <<- c(warn_bucket, conditionMessage(w))
        invokeRestart("muffleWarning")
      },
      message = function(m) {
        msg_bucket <<- c(msg_bucket, conditionMessage(m))
        invokeRestart("muffleMessage")
      }
    )
  }

  message("Running powerSim() for test term: ", test_term)
  ps <- with_warning_capture(
    powerSim(model_obj, test = fixed(test_term), nsim = nsim, alpha = ALPHA, progress = TRUE)
  )
  message("Running powerCurve() across subject breaks: ", paste(breaks, collapse = ", "))
  pc <- with_warning_capture(
    powerCurve(model_obj, test = fixed(test_term), along = along, breaks = breaks, nsim = nsim, alpha = ALPHA, progress = TRUE)
  )
  diag <- extract_warning_diagnostics(warn_bucket)

  list(
    power_sim = ps,
    power_curve = pc,
    diagnostics = diag,
    warnings = warn_bucket,
    messages = msg_bucket
  )
}

run_scenario <- function(beta_quadratic, random_slope_sd) {
  dat <- make_design(N_SUBJECT_BASE, TRIALS_PER_CONDITION)

  ## Simulate an inverted-U by combining linear and quadratic components.
  x_centered <- dat$contrast_num_c
  subject_intercepts <- rnorm(nlevels(dat$Subject), 0, 0.45)
  subject_slopes <- rnorm(nlevels(dat$Subject), 0, random_slope_sd)

  dat$gamma_power <- 2.4 +
    subject_intercepts[dat$Subject] +
    subject_slopes[dat$Subject] * x_centered +
    0.10 * x_centered +
    beta_quadratic * (x_centered^2) +
    rnorm(nrow(dat), 0, 1.20)

  ## Pilot-based primary model: random intercept.
  ## A random-slope sensitivity model is attempted; intercept-only is retained
  ## if the slope model is singular or fails.
  primary_model <- suppressWarnings(
    lmer(gamma_power ~ contrast_num_c + contrast_num_c2 + (1 | Subject), data = dat, REML = FALSE)
  )
  slope_model <- suppressWarnings(
    lmer(gamma_power ~ contrast_num_c + contrast_num_c2 + (1 + contrast_num_c | Subject), data = dat, REML = FALSE)
  )
  use_fallback <- isSingular(slope_model, tol = 1e-4)
  active_model <- primary_model
  if (!use_fallback) {
    active_model <- slope_model
  }

  ## Explicit inferential target aligns with H3 curvature term.
  test_term <- "contrast_num_c2"

  extended_model <- extend(active_model, along = "Subject", n = max(SUBJECT_BREAKS))
  out <- run_power_with_diagnostics(
    model_obj = extended_model,
    test_term = test_term,
    nsim = NSIM_ACTIVE,
    along = "Subject",
    breaks = SUBJECT_BREAKS
  )

  list(
    model = active_model,
    used_fallback = use_fallback,
    test_term = test_term,
    power_sim = out$power_sim,
    power_curve = out$power_curve,
    diagnostics = out$diagnostics
  )
}

scenario_grid <- expand.grid(
  effect_idx = seq_len(nrow(effect_scenarios)),
  slope_idx = seq_along(random_slope_sd_scenarios)
)

results <- vector("list", nrow(scenario_grid))
summary_rows <- vector("list", nrow(scenario_grid))
curve_rows <- list()

for (i in seq_len(nrow(scenario_grid))) {
  e_idx <- scenario_grid$effect_idx[i]
  s_idx <- scenario_grid$slope_idx[i]

  effect_label <- effect_scenarios$effect_label[e_idx]
  beta_quadratic <- effect_scenarios$beta_quadratic[e_idx]
  slope_sd <- random_slope_sd_scenarios[s_idx]
  message(
    "Scenario ", i, "/", nrow(scenario_grid),
    " | effect=", effect_label,
    " (beta_quadratic=", beta_quadratic, ")",
    " | random_slope_sd=", slope_sd
  )

  sim_out <- run_scenario(beta_quadratic = beta_quadratic, random_slope_sd = slope_sd)
  results[[i]] <- sim_out

  ps <- summary(sim_out$power_sim)
  summary_rows[[i]] <- data.frame(
    outcome = "gamma_power",
    effect_label = effect_label,
    beta_quadratic = beta_quadratic,
    random_slope_sd = slope_sd,
    test_term = sim_out$test_term,
    fallback_random_intercept = sim_out$used_fallback,
    power = ps$mean,
    lower = ps$lower,
    upper = ps$upper,
    mc_se = sqrt(ps$mean * (1 - ps$mean) / NSIM_ACTIVE),
    nsim = NSIM_ACTIVE,
    diagnostics_warnings = sim_out$diagnostics$n_warnings,
    diagnostics_convergence_rate = sim_out$diagnostics$convergence_rate,
    diagnostics_singular_rate = sim_out$diagnostics$singular_rate
  )

  pc <- summary(sim_out$power_curve)
  curve_rows[[i]] <- data.frame(
    outcome = "gamma_power",
    effect_label = effect_label,
    beta_quadratic = beta_quadratic,
    random_slope_sd = slope_sd,
    Subjects = pc$nlevels,
    Power = pc$mean,
    Lower = pc$lower,
    Upper = pc$upper,
    nsim = NSIM_ACTIVE
  )
}

scenario_summary <- do.call(rbind, summary_rows)
power_curve_df <- do.call(rbind, curve_rows)

print(scenario_summary)

## Save machine-readable outputs for manuscript supplement.
resolve_output_dir <- function() {
  preferred <- "/Volumes/g_psyplafor_methlab$/Students/Arne/GCP/figures/power_analysis"
  fallback <- file.path(getwd(), "_power_analysis", "outputs")
  if (dir.exists(dirname(preferred))) {
    return(preferred)
  }
  fallback
}
output_dir <- resolve_output_dir()
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
write.csv(
  scenario_summary,
  file = file.path(output_dir, "GCP_power_analysis_gamma_power_summary.csv"),
  row.names = FALSE
)
write.csv(
  power_curve_df,
  file = file.path(output_dir, "GCP_power_analysis_gamma_power_curve.csv"),
  row.names = FALSE
)
saveRDS(
  results,
  file = file.path(output_dir, "GCP_power_analysis_gamma_power_results.rds")
)

## Plot power curves across sensitivity scenarios.
png(
  file = file.path(output_dir, "GCP_power_analysis_gamma_power.png"),
  width = 2000, height = 1400, res = 300
)
ggplot(
  power_curve_df,
  aes(x = Subjects, y = Power, color = effect_label, linetype = factor(random_slope_sd))
) +
  geom_line() +
  geom_point(size = 1.8) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 1.5, alpha = 0.25) +
  geom_hline(yintercept = 0.90, linetype = "dashed", color = "grey30") +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent_format()) +
  labs(
    x = "Subjects",
    y = "Power",
    color = "Effect scenario",
    linetype = "Random-slope SD",
    title = "SIMR Power Curves: Gamma Power (contrast_num_c2)"
  ) +
  theme_minimal(base_size = 14)
dev.off()
