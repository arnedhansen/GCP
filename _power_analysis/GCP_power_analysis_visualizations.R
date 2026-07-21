## SIMR sensitivity analyses and custom power visualizations.
## The literature effect is fixed while pilot residual SD is varied.

script_dir <- {
  cmd <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd, value = TRUE)
  if (length(file_arg) > 0L) {
    dirname(sub("^--file=", "", file_arg[1]))
  } else {
    file.path(getwd(), "_power_analysis")
  }
}
source(file.path(script_dir, "GCP_power_analysis.R"))
# nolint start: object_usage_linter

ensure_packages(c("ggplot2", "scales"))
suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
})

run_residual_sensitivity <- function(
    label,
    outcome,
    dz,
    output_prefix,
    nsim = 1000L,
    subject_breaks = seq(10, 100, by = 10),
    seed = 123L,
    plot_only = FALSE) {
  gcp_root <- resolve_gcp_root()
  figure_dir <- file.path(gcp_root, "figures", "power_analysis")
  data_dir <- file.path(gcp_root, "data", "power_analysis")
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

  parameters <- load_variance_parameters(gcp_root, outcome, model_scope = "linear")
  scenario_labels <- c("Low", "Median", "High")
  scenario_sigma <- unname(parameters$residual_sd[c("low", "median", "high")])
  reference_sigma <- unname(parameters$residual_sd["median"])
  beta <- dz_to_linear_beta(dz, reference_sigma)
  csv_path <- file.path(data_dir, paste0(output_prefix, "_curve.csv"))

  if (plot_only) {
    if (!file.exists(csv_path)) {
      stop("Cached sensitivity curve not found: ", csv_path)
    }
    power_df <- utils::read.csv(csv_path, stringsAsFactors = FALSE)
  } else {
    rows <- lapply(seq_along(scenario_labels), function(i) {
      message(sprintf(
        "%s | %s residual SD = %.4f | d_z = %.3f | fixed beta = %.4f",
        label, scenario_labels[i], scenario_sigma[i], dz, beta
      ))
      model <- make_simr_model(
        n_subjects = max(subject_breaks),
        outcome_mean = parameters$outcome_mean,
        beta = beta,
        random_intercept_sd = unname(parameters$random_intercept_sd["median"]),
        residual_sd = scenario_sigma[i]
      )
      set.seed(seed + i - 1L)
      curve <- simr::powerCurve(
        model,
        test = simr::fixed("contrast_num_c", method = "t"),
        along = "Subject",
        breaks = subject_breaks,
        nsim = nsim,
        progress = TRUE
      )
      out <- extract_power_curve(curve, subject_breaks)
      out$scenario <- scenario_labels[i]
      out$residual_sd <- scenario_sigma[i]
      out
    })
    power_df <- do.call(rbind, rows)
    power_df$outcome <- label
    power_df$cohens_dz <- dz
    power_df$linear_beta <- beta
    power_df$reference_residual_sd <- reference_sigma
    utils::write.csv(power_df, csv_path, row.names = FALSE)
  }

  power_df$scenario <- factor(
    power_df$scenario,
    levels = scenario_labels,
    ordered = TRUE
  )
  legend_labels <- stats::setNames(
    sprintf(
      "%s residual SD: %.3f",
      scenario_labels,
      scenario_sigma
    ),
    scenario_labels
  )
  colors <- c(Low = "#1B9E77", Median = "#D95F02", High = "#7570B3")

  curve_plot <- ggplot2::ggplot(
    power_df,
    ggplot2::aes(
      x = .data$n_subjects,
      y = .data$power,
      color = .data$scenario,
      group = .data$scenario
    )
  ) +
    ggplot2::geom_line(linewidth = 0.9, linetype = "dotted") +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
      linewidth = 0.6,
      width = 1.6,
      alpha = 0.70
    ) +
    ggplot2::geom_point(size = 2.8) +
    ggplot2::geom_hline(
      yintercept = 0.90,
      linetype = "dashed",
      color = "grey45",
      linewidth = 0.7
    ) +
    ggplot2::scale_color_manual(
      values = colors,
      breaks = scenario_labels,
      labels = legend_labels
    ) +
    ggplot2::scale_x_continuous(breaks = subject_breaks) +
    ggplot2::scale_y_continuous(
      labels = scales::percent_format(accuracy = 1),
      limits = c(0, 1),
      breaks = c(0, 0.25, 0.50, 0.75, 0.90, 1)
    ) +
    ggplot2::labs(
      x = "Subjects",
      y = "Power",
      color = "Pilot uncertainty"
    ) +
    ggplot2::theme_classic(base_size = 18, base_family = "Arial") +
    ggplot2::theme(
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

  heatmap_plot <- ggplot2::ggplot(
    power_df,
    ggplot2::aes(
      x = factor(.data$n_subjects),
      y = .data$scenario,
      fill = .data$power
    )
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 1.1) +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.2f", .data$power)),
      color = "white",
      size = 3.52,
      fontface = "bold",
      family = "Arial"
    ) +
    ggplot2::scale_fill_gradientn(
      colours = c("#8E0F1F", "#F46D43", "#FEE08B", "#66BD63", "#0B6E3A"),
      values = c(0.00, 0.60, 0.79, 0.80, 1.00),
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.2),
      labels = sprintf("%.2f", seq(0, 1, by = 0.2))
    ) +
    ggplot2::scale_y_discrete(labels = legend_labels) +
    ggplot2::coord_fixed() +
    ggplot2::labs(x = "Subjects", y = NULL, fill = "Power") +
    ggplot2::theme_minimal(base_size = 16, base_family = "Arial") +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.title = ggplot2::element_text(size = 13),
      axis.text = ggplot2::element_text(size = 10),
      legend.position = "right",
      legend.title = ggplot2::element_text(size = 11),
      legend.text = ggplot2::element_text(size = 10)
    )

  grDevices::png(
    file.path(figure_dir, paste0(output_prefix, ".png")),
    width = 2200,
    height = 1400,
    res = 300
  )
  print(curve_plot)
  grDevices::dev.off()

  grDevices::png(
    file.path(figure_dir, paste0(output_prefix, "_heatmap.png")),
    width = 2500,
    height = 1400,
    res = 300
  )
  print(heatmap_plot)
  grDevices::dev.off()

  message("Saved custom visualizations for ", label, " to: ", figure_dir)
  invisible(power_df)
}

run_power_visualizations <- function(nsim = 1000L, plot_only = FALSE) {
  frequency <- run_residual_sensitivity(
    label = "Gamma peak frequency",
    outcome = "PeakFrequency",
    dz = 0.45,
    output_prefix = "GCP_power_analysis_frequency",
    nsim = nsim,
    seed = 223L,
    plot_only = plot_only
  )
  power <- run_residual_sensitivity(
    label = "Gamma peak power",
    outcome = "PeakAmplitude",
    dz = 0.52,
    output_prefix = "GCP_power_analysis_power",
    nsim = nsim,
    seed = 323L,
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
  run_power_visualizations(
    nsim = nsim,
    plot_only = "--plot-only" %in% args
  )
}
# nolint end
