## GCP GED Gamma Metrics — LME + Pilot Visualisation
##
## Part 1: LME computation (kept for future full-sample analysis)
## Part 2: Pilot-appropriate descriptive figures
##
## Dependencies: lme4, lmerTest, effectsize, ggplot2, dplyr, tidyr

library(lme4)
library(lmerTest)
library(effectsize)
library(ggplot2)
library(dplyr)
library(tidyr)

## ── Load data ────────────────────────────────────────────────────────
if (.Platform$OS.type == "windows") {
  csv_path <- "W:/Students/Arne/GCP/data/features/GCP_eeg_GED_gamma_metrics_trials.csv"
  fig_dir  <- "W:/Students/Arne/GCP/figures/stats/ged_gamma_metrics"
} else {
  csv_path <- "/Volumes/g_psyplafor_methlab$/Students/Arne/GCP/data/features/GCP_eeg_GED_gamma_metrics_trials.csv"
  fig_dir  <- "/Volumes/g_psyplafor_methlab$/Students/Arne/GCP/figures/stats/ged_gamma_metrics"
}
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

df <- read.csv(csv_path)

## Exclude S607 (extreme outlier: lambda1 = 5.15, dominates power variance)
exclude_subjects <- c(607)
df <- df[!df$Subject %in% exclude_subjects, ]
cat(sprintf("Excluded subjects: %s\n", paste(exclude_subjects, collapse = ", ")))

df$Subject   <- factor(df$Subject)
df$Condition <- factor(df$Condition, labels = c("25%", "50%", "75%", "100%"))

n_subj <- length(unique(df$Subject))
cat(sprintf("Loaded %d trials from %d subjects\n", nrow(df), n_subj))
cat(sprintf("Trials per condition: %s\n",
            paste(table(df$Condition), collapse = " / ")))

## ── Define metrics ───────────────────────────────────────────────────
metrics <- c("PeakFrequency", "PeakAmplitude", "BroadbandPower",
             "LowGammaPower", "HighGammaPower", "AUC", "FWHM",
             "Prominence", "SpectralCentroid", "PeakCount", "LoHiRatio")

metric_labels <- c("Peak Frequency [Hz]", "Peak Amplitude", "Broadband Power",
                    "Low-Gamma Power", "High-Gamma Power", "AUC Above Zero",
                    "FWHM [Hz]", "Peak Prominence", "Spectral Centroid [Hz]",
                    "Peak Count", "Low/High Ratio")

metrics_available <- metrics[metrics %in% names(df)]
labels_available  <- metric_labels[metrics %in% names(df)]

cond_colors <- c("25%" = "#FFD700", "50%" = "#FF8C00",
                 "75%" = "#DC143C", "100%" = "#2F2F2F")

## ======================================================================
##  PART 1: LME COMPUTATION (for future full-sample use)
## ======================================================================

## ── Model 1: Linear contrast effect ─────────────────────────────────
cat("\n========================================\n")
cat("MODEL 1: Linear contrast trend\n")
cat("  metric ~ Contrast + (1 | Subject)\n")
cat("========================================\n\n")

results_linear <- data.frame(
  Metric = character(), Beta = numeric(), SE = numeric(),
  t = numeric(), df = numeric(), p = numeric(),
  CI_lo = numeric(), CI_hi = numeric(), d = numeric(),
  stringsAsFactors = FALSE
)

for (i in seq_along(metrics_available)) {
  m <- metrics_available[i]
  d_tmp <- df[!is.na(df[[m]]), ]
  if (nrow(d_tmp) < 50) next

  d_tmp$Contrast_sc <- d_tmp$Contrast / 100
  fit <- lmer(reformulate("Contrast_sc + (1 | Subject)", response = m), data = d_tmp)
  s   <- summary(fit)
  cc  <- coef(s)["Contrast_sc", ]
  ci  <- confint(fit, parm = "Contrast_sc", method = "Wald")
  cohens_d <- cc["Estimate"] / sigma(fit)

  results_linear <- rbind(results_linear, data.frame(
    Metric = labels_available[i],
    Beta   = round(cc["Estimate"], 4),
    SE     = round(cc["Std. Error"], 4),
    t      = round(cc["t value"], 3),
    df     = round(cc["df"], 1),
    p      = cc["Pr(>|t|)"],
    CI_lo  = round(ci[1], 4),
    CI_hi  = round(ci[2], 4),
    d      = round(cohens_d, 3)
  ))

  cat(sprintf("%-22s  beta = %8.4f  SE = %7.4f  t(%5.1f) = %6.3f  p = %.4f  d = %6.3f\n",
              labels_available[i], cc["Estimate"], cc["Std. Error"],
              cc["df"], cc["t value"], cc["Pr(>|t|)"], cohens_d))
}

results_linear$p_adj <- p.adjust(results_linear$p, method = "holm")

cat("\n-- Holm-corrected p-values --\n")
for (r in 1:nrow(results_linear)) {
  cat(sprintf("%-22s  p_adj = %.4f\n",
              results_linear$Metric[r], results_linear$p_adj[r]))
}

## ── Model 2: Categorical condition (omnibus F-test) ──────────────────
cat("\n========================================\n")
cat("MODEL 2: Categorical condition (omnibus F)\n")
cat("  metric ~ Condition + (1 | Subject)\n")
cat("========================================\n\n")

results_categorical <- data.frame(
  Metric = character(), F_val = numeric(), num_df = numeric(),
  den_df = numeric(), p = numeric(), stringsAsFactors = FALSE
)

for (i in seq_along(metrics_available)) {
  m <- metrics_available[i]
  d_tmp <- df[!is.na(df[[m]]), ]
  if (nrow(d_tmp) < 50) next

  fit <- lmer(reformulate("Condition + (1 | Subject)", response = m), data = d_tmp)
  a   <- anova(fit, type = 3, ddf = "Satterthwaite")

  results_categorical <- rbind(results_categorical, data.frame(
    Metric = labels_available[i],
    F_val  = round(a["Condition", "F value"], 3),
    num_df = a["Condition", "NumDF"],
    den_df = round(a["Condition", "DenDF"], 1),
    p      = a["Condition", "Pr(>F)"]
  ))

  cat(sprintf("%-22s  F(%d, %5.1f) = %7.3f  p = %.4f\n",
              labels_available[i], a["Condition", "NumDF"],
              a["Condition", "DenDF"], a["Condition", "F value"],
              a["Condition", "Pr(>F)"]))
}

results_categorical$p_adj <- p.adjust(results_categorical$p, method = "holm")

cat("\n-- Holm-corrected p-values --\n")
for (r in 1:nrow(results_categorical)) {
  cat(sprintf("%-22s  p_adj = %.4f\n",
              results_categorical$Metric[r], results_categorical$p_adj[r]))
}

## ── Post-hoc pairwise comparisons for significant metrics ────────────
if (requireNamespace("emmeans", quietly = TRUE)) {
  library(emmeans)
  sig_metrics <- results_categorical$Metric[results_categorical$p_adj < 0.05]
  if (length(sig_metrics) > 0) {
    cat("\n========================================\n")
    cat("POST-HOC: Pairwise contrasts (significant metrics)\n")
    cat("========================================\n\n")
    for (i in seq_along(metrics_available)) {
      if (!labels_available[i] %in% sig_metrics) next
      m <- metrics_available[i]
      d_tmp <- df[!is.na(df[[m]]), ]
      fit <- lmer(reformulate("Condition + (1 | Subject)", response = m), data = d_tmp)
      emm <- emmeans(fit, "Condition")
      pw  <- pairs(emm, adjust = "holm")
      cat(sprintf("\n-- %s --\n", labels_available[i]))
      print(summary(pw))
    }
  }
}

## ── Save results tables ──────────────────────────────────────────────
write.csv(results_linear, file.path(fig_dir, "GCP_ged_gamma_LME_linear.csv"),
          row.names = FALSE)
write.csv(results_categorical, file.path(fig_dir, "GCP_ged_gamma_LME_categorical.csv"),
          row.names = FALSE)

## ======================================================================
##  PART 2: PILOT DESCRIPTIVE FIGURES
## ======================================================================

## ── Subject-level means ──────────────────────────────────────────────
df_subj <- df %>%
  group_by(Subject, Condition, Contrast) %>%
  summarise(across(all_of(metrics_available), ~mean(.x, na.rm = TRUE)),
            .groups = "drop")

## ── Figure 1: Pilot CRF Overview ────────────────────────────────────
## Multi-panel spaghetti plot: individual subjects + group mean +/- SEM
crf_metrics <- c("PeakFrequency", "SpectralCentroid", "PeakAmplitude",
                 "BroadbandPower", "LowGammaPower", "HighGammaPower",
                 "LoHiRatio", "AUC")
crf_labels  <- c("Peak Frequency [Hz]", "Spectral Centroid [Hz]",
                 "Peak Amplitude", "Broadband Power",
                 "Low-Gamma Power", "High-Gamma Power",
                 "Low/High Ratio", "AUC Above Zero")

crf_present <- crf_metrics %in% names(df_subj)
crf_metrics <- crf_metrics[crf_present]
crf_labels  <- crf_labels[crf_present]

df_crf <- df_subj %>%
  select(Subject, Condition, Contrast, all_of(crf_metrics)) %>%
  pivot_longer(cols = all_of(crf_metrics), names_to = "Metric", values_to = "Value") %>%
  filter(!is.na(Value))
df_crf$Metric <- factor(df_crf$Metric, levels = crf_metrics, labels = crf_labels)

df_crf_group <- df_crf %>%
  group_by(Metric, Condition, Contrast) %>%
  summarise(mu = mean(Value, na.rm = TRUE),
            sem = sd(Value, na.rm = TRUE) / sqrt(n()),
            .groups = "drop")

p_crf <- ggplot() +
  geom_line(data = df_crf,
            aes(x = Contrast, y = Value, group = Subject),
            color = "grey75", linewidth = 0.5, alpha = 0.7) +
  geom_point(data = df_crf,
             aes(x = Contrast, y = Value, color = Condition),
             size = 1.5, alpha = 0.6) +
  geom_ribbon(data = df_crf_group,
              aes(x = Contrast, ymin = mu - sem, ymax = mu + sem),
              fill = "black", alpha = 0.15) +
  geom_line(data = df_crf_group,
            aes(x = Contrast, y = mu),
            color = "black", linewidth = 1.4) +
  geom_point(data = df_crf_group,
             aes(x = Contrast, y = mu),
             color = "black", size = 3) +
  facet_wrap(~Metric, scales = "free_y", ncol = 4) +
  scale_color_manual(values = cond_colors) +
  scale_x_continuous(breaks = c(25, 50, 75, 100)) +
  labs(x = "Contrast [%]", y = NULL,
       title = "Pilot Contrast Response Functions (N = 9)",
       subtitle = "Grey = individual subjects | Black = group mean +/- SEM") +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "none",
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    strip.text = element_text(face = "bold", size = 11),
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 11, color = "grey40")
  )

ggsave(file.path(fig_dir, "GCP_pilot_CRF_overview.png"),
       p_crf, width = 14, height = 7, dpi = 300)

## ── Figure 2: Peak Frequency Detail ─────────────────────────────────
## Paired dot-plot with subject lines + violin
if ("PeakFrequency" %in% names(df_subj)) {
  freq_metric <- "PeakFrequency"
  freq_label  <- "Peak Frequency [Hz]"
} else {
  freq_metric <- "SpectralCentroid"
  freq_label  <- "Spectral Centroid [Hz]"
}

p_freq <- ggplot(df_subj, aes(x = Condition, y = .data[[freq_metric]])) +
  geom_violin(aes(fill = Condition), alpha = 0.25, color = NA, width = 0.8) +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", color = "grey40") +
  geom_line(aes(group = Subject), color = "grey60", linewidth = 0.6, alpha = 0.7) +
  geom_point(aes(color = Condition), size = 3, alpha = 0.85) +
  scale_fill_manual(values = cond_colors) +
  scale_color_manual(values = cond_colors) +
  labs(x = "Contrast", y = freq_label,
       title = paste0(freq_label, " by Contrast Level"),
       subtitle = sprintf("N = %d | Lines connect individual subjects", n_subj)) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    plot.title = element_text(size = 15, face = "bold"),
    plot.subtitle = element_text(size = 11, color = "grey40")
  )

ggsave(file.path(fig_dir, "GCP_pilot_peak_frequency.png"),
       p_freq, width = 7, height = 6, dpi = 300)

## ── Figure 3: Low/High Ratio Detail ─────────────────────────────────
if ("LoHiRatio" %in% names(df_subj)) {
  p_lohi <- ggplot(df_subj, aes(x = Condition, y = LoHiRatio)) +
    geom_violin(aes(fill = Condition), alpha = 0.25, color = NA, width = 0.8) +
    geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", color = "grey40") +
    geom_line(aes(group = Subject), color = "grey60", linewidth = 0.6, alpha = 0.7) +
    geom_point(aes(color = Condition), size = 3, alpha = 0.85) +
    scale_fill_manual(values = cond_colors) +
    scale_color_manual(values = cond_colors) +
    labs(x = "Contrast", y = "Low/High Gamma Ratio",
         title = "Low/High Gamma Ratio by Contrast Level",
         subtitle = sprintf("N = %d | Fraction of peaks in low gamma (30-50 Hz)", n_subj)) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "none",
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
      plot.title = element_text(size = 15, face = "bold"),
      plot.subtitle = element_text(size = 11, color = "grey40")
    )

  ggsave(file.path(fig_dir, "GCP_pilot_lohi_ratio.png"),
         p_lohi, width = 7, height = 6, dpi = 300)
}

## ── Figure 4: Effect Size Forest Plot ────────────────────────────────
## Clean forest plot showing Cohen's d for each metric (no significance stars)
results_linear$d_ci_lo <- results_linear$CI_lo / results_linear$SE * results_linear$d
results_linear$d_ci_hi <- results_linear$CI_hi / results_linear$SE * results_linear$d

p_forest <- ggplot(results_linear, aes(x = d, y = reorder(Metric, d))) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_errorbarh(aes(xmin = d_ci_lo, xmax = d_ci_hi),
                 height = 0.25, color = "grey40", linewidth = 0.6) +
  geom_point(size = 3.5, color = "black") +
  labs(x = "Cohen's d (linear contrast effect)", y = NULL,
       title = "Effect Sizes: Linear Contrast on Gamma Metrics",
       subtitle = sprintf("Pilot (N = %d) | LME: metric ~ Contrast + (1|Subject)", n_subj)) +
  theme_minimal(base_size = 14) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 12),
    plot.title = element_text(size = 15, face = "bold"),
    plot.subtitle = element_text(size = 11, color = "grey40")
  )

ggsave(file.path(fig_dir, "GCP_pilot_effect_sizes.png"),
       p_forest, width = 10, height = 6, dpi = 300)

## ── Figure 5: Key Power Metrics — Violin + Subject Lines ────────────
power_metrics <- c("PeakAmplitude", "BroadbandPower", "LowGammaPower", "HighGammaPower")
power_labels  <- c("Peak Amplitude", "Broadband Power", "Low-Gamma Power", "High-Gamma Power")
power_present <- power_metrics %in% names(df_subj)
power_metrics <- power_metrics[power_present]
power_labels  <- power_labels[power_present]

if (length(power_metrics) > 0) {
  df_pow <- df_subj %>%
    select(Subject, Condition, all_of(power_metrics)) %>%
    pivot_longer(cols = all_of(power_metrics), names_to = "Metric", values_to = "Value") %>%
    filter(!is.na(Value))
  df_pow$Metric <- factor(df_pow$Metric, levels = power_metrics, labels = power_labels)

  p_pow <- ggplot(df_pow, aes(x = Condition, y = Value)) +
    geom_violin(aes(fill = Condition), alpha = 0.2, color = NA, width = 0.8) +
    geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", color = "grey40") +
    geom_line(aes(group = Subject), color = "grey60", linewidth = 0.5, alpha = 0.7) +
    geom_point(aes(color = Condition), size = 2, alpha = 0.85) +
    facet_wrap(~Metric, scales = "free_y", ncol = 2) +
    scale_fill_manual(values = cond_colors) +
    scale_color_manual(values = cond_colors) +
    labs(x = "Contrast", y = "Metric Value",
         title = "Power Metrics by Contrast Level",
         subtitle = sprintf("N = %d | Lines connect individual subjects", n_subj)) +
    theme_minimal(base_size = 13) +
    theme(
      legend.position = "none",
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
      strip.text = element_text(face = "bold", size = 12),
      plot.title = element_text(size = 15, face = "bold"),
      plot.subtitle = element_text(size = 11, color = "grey40")
    )

  ggsave(file.path(fig_dir, "GCP_pilot_power_metrics.png"),
         p_pow, width = 10, height = 8, dpi = 300)
}

cat("\n========================================\n")
cat("Done. Results and figures saved to:\n")
cat(sprintf("  %s\n", fig_dir))
cat("========================================\n")
