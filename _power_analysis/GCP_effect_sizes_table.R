# GCP Effect Sizes Table
# Literature overview of grating contrast effects on gamma power and frequency.
# Effect sizes are derived within-subject dz, not values reported by the papers.
#
# Key output (data/power_analysis/):
#   GCP_effect_sizes_table.docx
#   GCP_effect_sizes_table.csv
#   GCP_effect_sizes_table_note.txt
#
# Run: Rscript GCP_effect_sizes_table.R

ensure_packages <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0L) {
    install.packages(missing, repos = "https://cloud.r-project.org")
  }
}

ensure_packages(c("flextable", "officer"))
suppressPackageStartupMessages({
  library(flextable)
  library(officer)
})

resolve_output_dirs <- function() {
  candidates <- c(
    "/Volumes/g_psyplafor_methlab$/Students/Arne/GCP/data/power_analysis",
    file.path(getwd(), "data", "power_analysis"),
    file.path(dirname(getwd()), "data", "power_analysis")
  )
  existing <- candidates[dir.exists(dirname(candidates)) | dir.exists(candidates)]
  if (length(existing) == 0L) {
    existing <- candidates[1]
  }
  unique(existing)
}

TABLE_FONT <- "Arial"
TABLE_FONT_SIZE <- 10
EM_DASH <- "\u2014"
DAGGER <- "\u2020"

# Derived dz = t / sqrt(N) = sqrt(F / N) from one-df within-subject tests.
# These papers do not report Cohen's d.
koch_dz <- sqrt(18.773 / 12)
murty_dz <- 4.20 / sqrt(12)
orekhova_power_dz <- sqrt(66.2 / 17)
orekhova_freq_dz <- sqrt(14.3 / 14)
# Karvat, Ofir, and Landau (2024): paired t(35) for 100% vs 50% contrast.
# t / sqrt(36) is 0.532 and 0.453; the simulation uses 0.52 and 0.45.
landau_power_dz <- 0.52
landau_freq_dz <- 0.45

fmt_dz <- function(x) {
  if (!is.finite(x)) {
    return(EM_DASH)
  }
  sprintf("%.2f", x)
}

make_row <- function(reference, n, nhp, conditions, stimulus, sf, dz) {
  n_txt <- if (is.na(n) || !nzchar(as.character(n))) {
    EM_DASH
  } else {
    as.character(n)
  }
  if (isTRUE(nhp) && n_txt != EM_DASH) {
    n_txt <- paste0(n_txt, DAGGER)
  }
  data.frame(
    Reference = reference,
    N = n_txt,
    `Contrast Condition [%]` = conditions,
    `Grating Stimulus` = stimulus,
    `Spatial Frequency (cpd)` = sf,
    `Effect Size` = fmt_dz(dz),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

gamma_power <- rbind(
  make_row("Bartoli et al. (2019)", 7, FALSE, "20, 50, 100", "Static sine-wave", "1", NA_real_),
  make_row("Hall et al. (2005)", 9, FALSE, "3, 6, 12, 24, 48, 96", "Static sine-wave", "3", NA_real_),
  make_row("Henrie & Shapley (2005)", 11, TRUE, "2, 3, 4, 6, 8, 13, 17, 23, 33, 47, 68, 96", "Drifting sine-waves", EM_DASH, NA_real_),
  make_row("Jia et al. (2013)", 5, TRUE, "1, 3.2, 6.4, 12.5, 25, 50, 100", "Drifting sinusoidals", "1", NA_real_),
  make_row("Koch et al. (2009)", 12, FALSE, "5, 9, 24, 54, 91", "Concentric dynamic", "1", koch_dz),
  make_row("Murty et al. (2018)", 12, FALSE, "0, 12.5, 25, 37.5, 50, 62.5, 75, 87.5, 100", "Full-screen sinusoidal", "0.5 to 8", murty_dz),
  make_row("Orekhova et al. (2020)", 17, FALSE, "50, 100", "Concentric dynamic sine-wave", "1.66", orekhova_power_dz),
  make_row("Perry et al. (2015)", 18, FALSE, "40.9 to 100 (continuous)", "Annular square-wave", "3", NA_real_),
  make_row("Schadow et al. (2007)", 21, FALSE, "5, 25, 50", "Static sine-wave", "5", NA_real_),
  make_row("Van Pelt et al. (2018)", 158, FALSE, "50, 100", "Concentric dynamic sine-wave", "3", NA_real_),
  make_row("Karvat et al. (2024)", 36, FALSE, "50, 100", "Annular square-wave", "3", landau_power_dz)
)

gamma_frequency <- rbind(
  make_row("Bartoli et al. (2019)", 7, FALSE, "20, 50, 100", "Static sine-wave", "1", NA_real_),
  make_row("Hadjipapas et al. (2015)", 9, FALSE, "20, 36, 48, 66, 96", "Static square-wave", "3", NA_real_),
  make_row("Jia et al. (2013)", 5, TRUE, "1, 3.2, 6.4, 12.5, 25, 50, 100", "Drifting sinusoidals", "1", NA_real_),
  make_row("Lowet et al. (2015)", 1, TRUE, "6.1, 9.7, 16.3, 35.9, 50.3, 72", "Static square-wave", "2", NA_real_),
  make_row("Orekhova et al. (2020)", 17, FALSE, "50, 100", "Concentric dynamic sine-wave", "1.66", orekhova_freq_dz),
  make_row("Perry et al. (2015)", 18, FALSE, "40.9 to 100 (continuous)", "Annular square-wave", "3", NA_real_),
  make_row("Ray & Maunsell (2010)", 2, TRUE, "0, 1.6, 3.1, 6.2, 12.5, 25, 50, 100", "Gabor patches", "4", NA_real_),
  make_row("Roberts et al. (2013)", 2, TRUE, "2.5, 3.7, 6.1, 9.7, 16.3, 35.9, 50.3, 72", "Static square-wave", "2", NA_real_),
  make_row("Van Pelt et al. (2018)", 158, FALSE, "50, 100", "Concentric dynamic sine-wave", "3", NA_real_),
  make_row("Karvat et al. (2024)", 36, FALSE, "50, 100", "Annular square-wave", "3", landau_freq_dz)
)

TABLE_NOTE <- paste(
  "Table 1. Overview of studies reporting effects of grating contrast on gamma power and frequency.",
  "Contrast values are those stated by each paper and are not assumed to be Michelson unless that metric is named in the source.",
  "Effect Size reports derived within-subject dz, not a statistic published as Cohen's d.",
  "For one-df tests, dz = t / sqrt(N) = sqrt(F / N): Koch F(1,11) = 18.773, N = 12, dz = 1.25;",
  "Murty human EEG slow-gamma contrast slope t(11) = 4.20, N = 12, dz = 1.21;",
  "Orekhova power F(1,16) = 66.2, N = 17, dz = 1.97;",
  "Orekhova frequency F(1,13) = 14.3, implying N = 14 for that test, dz = 1.01;",
  "Karvat et al. (2024) paired t(35) = 3.19 for gamma power and t(35) = 2.72 for gamma peak cycle, 100% versus 50% contrast, entered as dz = 0.52 and dz = 0.45.",
  "Among convertible one-df derived dz values, these are the smallest in each section and are the planning inputs.",
  "Omnibus tests with more than one numerator df (Bartoli; Schadow) are not converted to a pairwise dz.",
  "Hadjipapas et al. (2015) is retained, but a previously listed dz of 2.14 is omitted because it used a between-subject conversion of a within-subject F.",
  paste0(DAGGER, " Nonhuman primate studies: N is the number of animals; inferential statistics typically rest on recording sites (Henrie 68 sites; Jia 90 sites; Ray 23 and 59 electrodes)."),
  "Lowet et al. (2015) reanalysed LFP spectra from Roberts et al. (2013) (one monkey illustrated; six contrast levels with a spectral peak).",
  "Hadjipapas et al. (2015) compared a new human MEG sample with the Roberts monkey dataset.",
  "Schadow et al. (2007) measured the early evoked gamma-band response rather than sustained induced gamma.",
  "Murty et al. (2018) recruited 19 participants; contrast analyses used N = 12, with spatial frequency selected per participant from 0.5 to 8 cpd (2 cpd is the macaque value).",
  "Koch et al. (2009) used concentric dynamic gratings whose luminance profile is not stated.",
  "Ray & Maunsell (2010) presented eight contrasts as Gabor patches; frequency analyses used 25, 50, and 100%.",
  "Henrie & Shapley (2005) reported Rayleigh contrast, numerically equivalent to Michelson for these gratings.",
  "Perry et al. (2015) mapped amplitude and frequency continuously from 40.9% to 100% Michelson contrast and is listed in both sections.",
  "Karvat et al. (2024) used an annular square-wave grating at 3 cpd; N = 36 is the EEG contrast comparison (df = 35)."
)

style_effect_size_ft <- function(ft, col_widths = NULL) {
  ft <- theme_booktabs(ft)
  rule_border <- fp_border(color = "#666666", width = 0.75)
  ft <- hline_top(ft, border = rule_border, part = "header")
  ft <- hline_bottom(ft, border = rule_border, part = "header")
  ft <- hline_bottom(ft, border = rule_border, part = "body")
  ft <- font(ft, fontname = TABLE_FONT, part = "all")
  ft <- fontsize(ft, size = TABLE_FONT_SIZE, part = "all")
  ft <- bold(ft, part = "header")
  ft <- padding(ft, padding = 3, part = "all")
  ft <- align(ft, align = "left", j = 1, part = "all")
  if (ncol(ft$body$dataset) > 1) {
    ft <- align(ft, align = "center", j = 2:ncol(ft$body$dataset), part = "all")
  }
  ft <- valign(ft, valign = "center", part = "all")
  ft <- set_table_properties(ft, layout = "fixed")
  if (!is.null(col_widths) && length(col_widths) == ncol(ft$body$dataset)) {
    ft <- width(ft, j = seq_along(col_widths), width = col_widths)
  } else {
    ft <- autofit(ft)
  }
  ft
}

make_contrast_flextable <- function(df) {
  ft <- flextable(df)
  style_effect_size_ft(
    ft,
    col_widths = c(1.8, 0.40, 1.55, 1.45, 0.70, 0.55)
  )
}

add_section <- function(doc, title, ft) {
  doc <- body_add_fpar(doc, fpar(
    ftext(title, prop = fp_text(bold = TRUE, font.family = TABLE_FONT, font.size = TABLE_FONT_SIZE))
  ))
  doc <- body_add_flextable(doc, ft)
  body_add_par(doc, "")
}

export_table <- function(out_dir, stacked) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  doc <- read_docx()
  doc <- add_section(
    doc,
    "Studies reporting on the effect of contrast on gamma power",
    make_contrast_flextable(gamma_power)
  )
  doc <- add_section(
    doc,
    "Studies reporting on the effect of contrast on gamma frequency",
    make_contrast_flextable(gamma_frequency)
  )
  doc <- body_add_fpar(doc, fpar(
    ftext(TABLE_NOTE, prop = fp_text(font.family = TABLE_FONT, font.size = 9, italic = TRUE))
  ))

  docx_path <- file.path(out_dir, "GCP_effect_sizes_table.docx")
  csv_path <- file.path(out_dir, "GCP_effect_sizes_table.csv")
  note_path <- file.path(out_dir, "GCP_effect_sizes_table_note.txt")
  print(doc, target = docx_path)
  utils::write.csv(stacked, csv_path, row.names = FALSE)
  writeLines(TABLE_NOTE, note_path)
  message("Saved: ", docx_path)
  message("Saved: ", csv_path)
  message("Saved: ", note_path)
}

stacked <- rbind(
  cbind(Section = "Gamma power", gamma_power, stringsAsFactors = FALSE),
  cbind(Section = "Gamma frequency", gamma_frequency, stringsAsFactors = FALSE)
)

for (out_dir in resolve_output_dirs()) {
  tryCatch(
    export_table(out_dir, stacked),
    error = function(e) message("Could not write to ", out_dir, ": ", conditionMessage(e))
  )
}

message("Derived dz used in the table:")
message(sprintf("  Koch power:     %.4f", koch_dz))
message(sprintf("  Murty power:    %.4f", murty_dz))
message(sprintf("  Orekhova power: %.4f", orekhova_power_dz))
message(sprintf("  Orekhova freq:  %.4f", orekhova_freq_dz))
message(sprintf("  Landau power:   %.4f", landau_power_dz))
message(sprintf("  Landau freq:    %.4f", landau_freq_dz))
